#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ParaSparse.h"

const ParaSparse PS_DEFAULT = {.na=0, .nd=0, .N=0, .Ne=0, .length=0, .i=NULL, .j=NULL, .Mij=NULL,
							.rank=0, .size=0, .comm=NULL, .send_counts=NULL, .recv_counts=NULL,
							.send_displs=NULL, .recv_displs=NULL, .N_send=0, .N_recv=0,
							.i_send=NULL, .i_recv=NULL, .y_send=NULL, .y_recv=NULL,
							.analyzed=0, .finalized=0};	

char grav1d_name[] = "grav1d";
char grav2d_name[] = "grav2d";
char sine_name[] = "sine";

void PS_generate_sine(ParaSparse *M, int *N, MPI_Comm comm)
{
	int i,j,n;
	
	*M = PS_DEFAULT;
	
	M->N = *N;
	M->comm = comm;
	MPI_Comm_size(M->comm, &(M->size));
	MPI_Comm_rank(M->comm, &(M->rank));
	
	M->na = (M->rank) * (*N/(M->size));
	if(M->rank == M->size-1)
		M->nd = *N - M->na;
	else
		M->nd = ((M->rank)+1)*(*N/(M->size)) - M->na;
	
	for(n = 0; n < M->nd; n++)
		PS_add_entry(M, n+M->na, n+M->na, 3.0+sin(n+M->na));
	
	for(i = 0; i < *N; i++)
		for(j = 0; j < *N; j++)
			if(j < i && i%(j+1) == 0 
			   && ((i>=M->na && i<M->na+M->nd) || (j>=M->na && j<M->na+M->nd)))
				PS_add_entry(M, i, j, 0.1*sin(i + (*N)*j));
}

void PS_generate_grav1d(ParaSparse *M, int *N, MPI_Comm comm)
{
	int i,j,n;
	
	*M = PS_DEFAULT;
	
	M->N = *N;
	M->comm = comm;
	MPI_Comm_size(M->comm, &(M->size));
	MPI_Comm_rank(M->comm, &(M->rank));
	
	M->na = (M->rank) * (*N/(M->size));
	if(M->rank == M->size-1)
		M->nd = *N - M->na;
	else
		M->nd = ((M->rank)+1)*(*N/(M->size)) - M->na;
	
	for(n = 0; n < M->nd; n++)
		PS_add_entry(M, n+M->na, n+M->na, -2.0);
	
	for(i = 0; i < *N; i++)
		for(j = 0; j < *N; j++)
			if(j == i-1 
			   && ((i>=M->na && i<M->na+M->nd) || (j>=M->na && j<M->na+M->nd)))
				PS_add_entry(M, i, j, 1.0);
	
}

void PS_generate_grav2d(ParaSparse *M, int *N, MPI_Comm comm)
{
	int i,j,n, Nx;
	int x1, y1, x2, y2;
	
	*M = PS_DEFAULT;
	
	Nx = (int) sqrt( *N);
	*N = Nx * Nx;
	M->N = *N;
	M->comm = comm;
	MPI_Comm_size(M->comm, &(M->size));
	MPI_Comm_rank(M->comm, &(M->rank));
	
	M->na = (M->rank) * (*N/(M->size));
	if(M->rank == M->size-1)
		M->nd = *N - M->na;
	else
		M->nd = ((M->rank)+1)*(*N/(M->size)) - M->na;
	
//	for(n = 0; n < M->nd; n++)
//		PS_add_entry(M, n+M->na, n+M->na, -4.0);
	
	for(i = 0; i < *N; i++)
		for(j = 0; j < *N; j++)
		{
			x1 = i % Nx;
			y1 = i / Nx;
			x2 = j % Nx;
			y2 = j / Nx;
			
			if(j <= i && ((i>=M->na && i<M->na+M->nd) || (j>=M->na && j<M->na+M->nd)))
			{
				if(i == j)
					PS_add_entry(M, i, j, -4.0);
				else if(((x1 == x2) && (y1 == y2+1 || y1 == y2-1))
						|| ((y1 == y2) && (x1 == x2+1 || x1 == x2-1)))
					PS_add_entry(M, i, j, 1.0);
			}
		}
	
}

void PS_add_entry(ParaSparse *M, int i, int j, double Mij)
{
	M->finalized = 0;
	M->analyzed = 0;
	
	if(i < j)
	{
		int temp = i;
		i = j;
		j = temp;
	}
	
	if(M->Ne > 0)
	{
		if(M->Ne >= M->length)
		{
			while(M->Ne >= M->length)
				M->length *= 2;
			
			M->i = (int *) realloc(M->i, M->length * sizeof(int));
			M->j = (int *) realloc(M->j, M->length * sizeof(int));
			M->Mij = (double *) realloc(M->Mij, M->length * sizeof(double));
		}
		
		M->i[M->Ne] = i;
		M->j[M->Ne] = j;
		M->Mij[M->Ne] = Mij;
		(M->Ne)++;
	}
	
	else
	{
		M->Ne = 1;
		M->length = 16;
		M->i = (int *) malloc((M->length) * sizeof(int));
		M->j = (int *) malloc((M->length) * sizeof(int));
		M->Mij = (double *) malloc((M->length) * sizeof(double));
		
		M->i[0] = i;
		M->j[0] = j;
		M->Mij[0] = Mij;
	}
}

void PS_finalize(ParaSparse *M)
{
	if(M->length > M->Ne)
	{
		M->length = M->Ne;
	
		M->i = (int *) realloc(M->i, M->length * sizeof(int));
		M->j = (int *) realloc(M->j, M->length * sizeof(int));
		M->Mij = (double *) realloc(M->Mij, M->length * sizeof(double));
	}
	
	M->finalized = 1;
}

void PS_analyze(ParaSparse *M)
{
	if(!(M->finalized))
		PS_finalize(M);
	
	int i,j,n,r,ind;
	
	int N = M->N;
	int *nas, *nds, *full_count;
	
	nas = (int *) malloc(M->size * sizeof(int));
	nds = (int *) malloc(M->size * sizeof(int));
	
	M->send_counts = (int *) malloc(M->size * sizeof(int));
	M->recv_counts = (int *) malloc(M->size * sizeof(int));
	M->send_displs = (int *) malloc(M->size * sizeof(int));
	M->recv_displs = (int *) malloc(M->size * sizeof(int));
	
	MPI_Allgather(&(M->na), 1, MPI_INT, nas, 1, MPI_INT, M->comm);
	MPI_Allgather(&(M->nd), 1, MPI_INT, nds, 1, MPI_INT, M->comm);
	
	for(n=0; n<M->size; n++)
	{	
		M->send_counts[n] = 0;
		M->recv_counts[n] = 0;
		M->send_displs[n] = 0;
		M->recv_displs[n] = 0;
	}
	
	full_count = (int *) malloc(N * sizeof(int));
	for(n=0; n<N; n++)
		full_count[n] = 0;
	
	for(n=0; n < M->Ne; n++)
	{
		i = M->i[n];
		j = M->j[n];
		
		if(i >= M->na && i < M->na + M->nd)
			full_count[j]++;
		if(j >= M->na && j < M->na + M->nd && j != i)
			full_count[i]++;
	}
	
	r = 0;
	M->N_send = 0;
	
//	printf("R%d: full_count.\n", M->rank);
//	for(n=0; n<N; n++)
//		printf("%d ", full_count[n]);
//	printf("\n");
	
	for(n=0; n<N; n++)
	{
		if(n >= nas[r] + nds[r])
			r++;
		
		if(full_count[n] > 0 && r != M->rank)
			(M->send_counts[r])++;
	}
	
	for(r = 0; r < M->size; r++)
	{
		M->N_send += M->send_counts[r];
		
		if(r > 0)
			M->send_displs[r] = M->send_displs[r-1] + M->send_counts[r-1];
	}
	
	M->i_send = (int *) malloc(M->N_send * sizeof(int));
	M->y_send = (double *) malloc(M->N_send * sizeof(double));
	
	for(n = 0; n < M->N_send; n++)
		M->i_send[n] = -1;
	
	for(n = 0; n < M->Ne; n++)
	{
		i = M->i[n];
		j = M->j[n];
		
		if(j < M->na)
		{
			r = j;
			j = i;
			i = r;
		}
		
		if(i < M->na || i >= M->na + M->nd)
		{
			for(r = 0; r < M->size; r++)
				if(i < nas[r] + nds[r])
					break;
			
			for(ind = 0; ind < M->send_counts[r]; ind++)
			{
				if(M->i_send[M->send_displs[r]+ind] == -1)
				{
					M->i_send[M->send_displs[r]+ind] = i;
					break;
				}
				if(M->i_send[M->send_displs[r]+ind] == i)
					break;
			}
		}
	}
	
	MPI_Alltoall(M->send_counts, 1, MPI_INT, M->recv_counts, 1, MPI_INT, M->comm);
	
	M->N_recv = 0;
	for(r = 0; r < M->size; r++)
	{
		M->N_recv += M->recv_counts[r];
		
		if(r > 0)
			M->recv_displs[r] = M->recv_displs[r-1] + M->recv_counts[r-1];
	}
	
	M->i_recv = (int *) malloc(M->N_recv * sizeof(int));
	M->y_recv = (double *) malloc(M->N_recv * sizeof(double));
	
	MPI_Alltoallv(M->i_send, M->send_counts, M->send_displs, MPI_INT,
				  M->i_recv, M->recv_counts, M->recv_displs, MPI_INT, M->comm);
	
	for(n=0; n < M->N_recv; n++)
	{
		if(M->i_recv[n] < M->na)
			printf("Rank %d received low index %d: na=%d nd=%d\n", M->rank, M->i_recv[n], M->na, M->nd);
		else if(M->i_recv[n] >= M->na + M->nd)
			printf("Rank %d received high index %d: na=%d nd=%d\n", M->rank, M->i_recv[n], M->na, M->nd);
	}
	
	free(nas);
	free(nds);
	free(full_count);
	
	M->analyzed = 1;
}

void PS_multiply(ParaSparse *M, double *x, double *y)
{
	if(!(M->analyzed))
		PS_analyze(M);
	
	int i,j,n,r;
	int na = M->na;
	int nd = M->nd;
	
	for(n=0; n<nd; n++)
		y[n] = 0;
	
	for(n=0; n<M->N_send; n++)
		M->y_send[n] = 0;
	
	for(n = 0; n < M->Ne; n++)
	{
		i = M->i[n];
		j = M->j[n];
		
		if(i >= na && i < na+nd && j >= na && j <= na+nd)
		{
			y[i-na] += M->Mij[n] * x[j-na];
			if(i != j)
				y[j-na] += M->Mij[n] * x[i-na];
		}
		
		else
		{
			if(j < na)
			{
				r = j;
				j = i;
				i = r;
			}
			
			//TODO: Make this not suck --> Lookup table?
			for(r = 0; r < M->N_send; r++)
				if(M->i_send[r] == i)
				{
					M->y_send[r] += M->Mij[n] * x[j-na];
					break;
				}
		}
	}
	
//	printf("R%d sending %d.\n", M->rank, M->N_send);
//	printf("R%d sending:", M->rank);
//	for(n=0; n<M->size; n++)
//		printf(" %d to %d,", M->send_counts[n], n);
//	printf("\n");
//	for(n=0; n<M->N_send; n++)
//	{	
//		printf("R%d sent: %d %lg\n", M->rank, M->i_send[n], M->y_send[n]);
//	}
	
	MPI_Alltoallv(M->y_send, M->send_counts, M->send_displs, MPI_DOUBLE,
				  M->y_recv, M->recv_counts, M->recv_displs, MPI_DOUBLE, M->comm);
	
//	printf("R%d receiving %d.\n", M->rank, M->N_recv);
//	printf("R%d receiving:", M->rank);
//	for(n=0; n<M->size; n++)
//		printf(" %d from %d,", M->recv_counts[n], n);
//	printf("\n");
	for(n=0; n<M->N_recv; n++)
	{	
		y[M->i_recv[n] - na] += M->y_recv[n];
//		printf("R%d received: %d %lg\n", M->rank, M->i_recv[n], M->y_recv[n]);
	}
}

void PS_printM(ParaSparse *M)
{
	int r,n;
	
	for(r = 0; r < M->size; r++)
	{
		if(r == M->rank)
		{	
			printf("Rank %d of %d: na=%d nd=%d\n", M->rank, M->size, M->na, M->nd);
			for(n = 0; n < M->Ne; n++)
				printf("%d %d %lg\n",M->i[n],M->j[n],M->Mij[n]);
		}
		MPI_Barrier(M->comm);
	}
	
	if(M->rank == 0)
		printf("\n");
}

void PS_free(ParaSparse *M)
{	
	free(M->i);
	free(M->j);
	free(M->Mij);
	free(M->send_counts);
	free(M->recv_counts);
	free(M->send_displs);
	free(M->recv_displs);
	free(M->i_send);
	free(M->i_recv);
	free(M->y_send);
	free(M->y_recv);
}

void PS_printV(double *v, int na, int nd, MPI_Comm comm)
{
	int r,n,rank,size;
	
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);
	
	for(r = 0; r < size; r++)
	{
		if(r == rank)
		{	
			printf("Rank %d of %d:\n", rank, size);
			for(n = 0; n < nd; n++)
				printf("%d %lg\n", n+na, v[n]);
		}
		MPI_Barrier(comm);
	}
	
	fflush(stdout);
	
	if(rank == 0)
		printf("\n");
}

double PS_dot(double *x, double *y, int nd, MPI_Comm comm)
{
	int i;
	double xy = 0.0;
	double xy_out;
	
	for(i=0; i<nd; i++)
		xy += x[i]*y[i];
	
	MPI_Allreduce(&xy, &xy_out, 1, MPI_DOUBLE, MPI_SUM, comm);
	
	return xy_out;
}

void PS_get_diag(ParaSparse *A, double *d)
{
	int n;
	
	int *nas, *nds;
	double *diag_loc;
	
	nas = (int *) malloc(A->size * sizeof(int));
	nds = (int *) malloc(A->size * sizeof(int));
	diag_loc = (double *) malloc(A->nd * sizeof(double));
	
	for(n = 0; n < A->nd; n++)
		diag_loc[n] = 0; //TODO: memset?
	
	for(n = 0; n < A->Ne; n++)
		if(A->i[n] == A->j[n])
			diag_loc[A->i[n] - A->na] += A->Mij[n];
	
	MPI_Allgather(&(A->na), 1, MPI_INT, nas, 1, MPI_INT, A->comm);
	MPI_Allgather(&(A->nd), 1, MPI_INT, nds, 1, MPI_INT, A->comm);
	
	MPI_Allgatherv(diag_loc, A->nd, MPI_DOUBLE,
				   d, nds, nas, MPI_DOUBLE, A->comm);
	
	free(nas);
	free(nds);
	free(diag_loc);
}

void PS_get_precond(ParaSparse *A, ParaSparse *C)
{
	int i,j,n;
	double *diag;
	
	*C = PS_DEFAULT;

	C->comm = A->comm;
	C->rank = A->rank;
	C->size = A->size;
	C->N = A->N;
	C->na = A->na;
	C->nd = A->nd;
	
	diag = (double *) malloc(A->N * sizeof(double));
	PS_get_diag(A, diag);
	
	for(n = 0; n < A->nd; n++)
		PS_add_entry(C, A->na + n, A->na + n, 1.0/diag[A->na + n]);
	
	for(n = 0; n < A->Ne; n++)
	{
		i = A->i[n];
		j = A->j[n];
		
		if(i != j)
			PS_add_entry(C, i, j, -(A->Mij[n]) / (diag[i]*diag[j]));
	}
	
	free(diag);
}

void PS_bcg_iter(ParaSparse *A, ParaSparse *C, double *r, double *p, double *x)
{
	int i;
	double r2, pAp, alpha, beta;
	double *rold, *pold, *Ap, *CAp;
	
	int nd = A->nd;
	
	//TODO: declare/allocate these in PS_bcg
	rold = (double *) malloc(nd * sizeof(double));
	pold = (double *) malloc(nd * sizeof(double));
	Ap = (double *) malloc(nd * sizeof(double));
	CAp = (double *) malloc(nd * sizeof(double));
	
	for(i=0; i<nd; i++)
	{
		rold[i] = r[i]; //TODO: memcpy?
		pold[i] = p[i]; //TODO: memcpy?
	}
	
	r2 = PS_dot(rold, rold, nd, A->comm);
	PS_multiply(A, pold, Ap);
	PS_multiply(C, Ap, CAp);
	pAp = PS_dot(pold, CAp, nd, A->comm);
	alpha = r2/pAp;
	
	for(i=0; i<nd; i++)
		r[i] = rold[i] - alpha * CAp[i];
	
	beta = PS_dot(r, r, nd, A->comm) / r2;
	
	for(i=0; i<nd; i++)
	{
		p[i] = r[i] + beta * pold[i];
		x[i] += alpha * pold[i];
	}
	
	free(rold);
	free(pold);
	free(Ap);
	free(CAp);
}

void PS_bcg(ParaSparse *A, double *b, double *x)
{
	int i;
	
	double res0, err;
	double *Ax, *CAx, *Cb, *r, *p;
	int nd = A->nd;
	
	Ax = (double *) malloc(nd * sizeof(double));
	CAx = (double *) malloc(nd * sizeof(double));
	Cb = (double *) malloc(nd * sizeof(double));
	r = (double *) malloc(nd * sizeof(double));
	p = (double *) malloc(nd * sizeof(double));
	
	ParaSparse C;
	PS_get_precond(A, &C);
	
	PS_multiply(&C, b, Cb);

	for(i=0; i<nd; i++)
		x[i] = Cb[i];
	
	PS_multiply(A, x, Ax);
	PS_multiply(&C, Ax, CAx);
	
	for(i=0; i<nd; i++)
	{
		r[i] = Cb[i] - CAx[i];
		p[i] = r[i];  //TODO: memcpy?
	}
	
	res0 = PS_dot(Cb, Cb, nd, A->comm);
	err = 1.0;
	
	i = 0;
	while(err > TOL && i < 2*(A->N))
	{
		PS_bcg_iter(A, &C, r, p, x);
//		PS_printV(x, A->na, nd, A->comm);
		err = sqrt(PS_dot(r,r,nd,A->comm) / res0);
		
		if(A->rank == 0)
			printf("err = %e iter = %d\n",err,i);
		
		i++;
	}
	
	PS_free(&C);
	
	free(Ax);
	free(CAx);
	free(Cb);
	free(r);
	free(p);
	
}
