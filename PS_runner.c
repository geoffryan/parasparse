#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "PS_runner.h"
#include "ParaSparse.h"

void run_PS(int argc, char *argv[])
{
	if(argc < 4)
	{
		printf("Please give problem, type (solve/multiply), and a matrix size.\n");
		return;
	}
	
	int i;
	int N = atoi(argv[3]);
	double *x, *y;
	clock_t ticks;
	
	if(N <= 0)
	{
		printf("The size must be positive numbnuts!\n");
		return;
	}
	
	ParaSparse A;
	
	if(strcmp(argv[1], sine_name) == 0)
		PS_generate_sine(&A, &N, MPI_COMM_WORLD);
	
	else if(strcmp(argv[1], grav1d_name) == 0)
		PS_generate_grav1d(&A, &N, MPI_COMM_WORLD);
	
	else if(strcmp(argv[1], grav2d_name) == 0)
		PS_generate_grav2d(&A, &N, MPI_COMM_WORLD);
	
	else if(strcmp(argv[1], grav3d_name) == 0)
		PS_generate_grav3d(&A, &N, MPI_COMM_WORLD);
	
	else
	{
		if(A.rank == 0) 
			printf("You call that a name numbnuts?\n");
		return;
	}
	
	if(A.N <= 16)
		PS_printM(&A);
	
	x = (double *) malloc(A.nd * sizeof(double));
	y = (double *) malloc(A.nd * sizeof(double));
	
	if(strcmp(argv[2], "multiply") == 0)
	{
		for(i=0; i<A.nd; i++)
			x[i] = A.na + i+1;
		
		if(A.rank == 0)
			printf("\nMultiplying a %dx%d %s system with %d processes.\n", N, N, argv[1], A.size);
		
		PS_analyze(&A);
		MPI_Barrier(MPI_COMM_WORLD);
		if(A.rank == 0)
			ticks = clock();
		
		PS_multiply(&A, x, y);
		
		MPI_Barrier(MPI_COMM_WORLD);
		if(A.rank == 0)
			ticks = clock() - ticks;
		
		double dot = PS_dot(y,y,A.nd,A.comm);
		
		if(A.rank == 0)
		{
			char filename[51];
			snprintf(filename, 51, "time_%s_mult.out", argv[1]);
			FILE* timefile = fopen(filename, "a");
			fprintf(timefile, "%d %d %lg %d %lg\n", N, A.size, dot, (int) ticks, ((double)ticks)/CLOCKS_PER_SEC);
			fclose(timefile);
		}
		
		if(A.rank == 0) 
			printf("<y,y> = %lg\n", dot);
	}
	
	else if(strcmp(argv[2], "solve") == 0)
	{
		if(A.rank == 0)
			printf("\nSetting up the problem, just be patient.\n");
		for(i=0; i<A.nd; i++)
			x[i] = A.na + i+1;
		PS_multiply(&A, x, y);
		for(i=0; i<A.nd; i++)
			x[i] = 1.0;
		
		if(A.rank == 0)
			printf("Solving a %dx%d %s system with %d processes.\n\n", N, N, argv[1], A.size);
	
		MPI_Barrier(MPI_COMM_WORLD);
		if(A.rank == 0)
			ticks = clock();
		
		int iters = PS_bcg(&A, y, x);
		
		MPI_Barrier(MPI_COMM_WORLD);
		if(A.rank == 0)
			ticks = clock() - ticks;
		
		double dot = PS_dot(x,x,A.nd,A.comm);
		
		if(A.rank == 0)
		{
			char filename[51];
			snprintf(filename, 51, "time_%s_solve.out", argv[1]);
			FILE* timefile = fopen(filename, "a");
			fprintf(timefile, "%d %d %d %lg %d %lg\n", N, A.size, iters, dot, (int) ticks, ((double)ticks)/CLOCKS_PER_SEC);
			fclose(timefile);
		}
		
		if(A.rank == 0) 
			printf("<x,x> = %lg\n", dot);
	}
	
	else
	{
		if(A.rank == 0) 
			printf("\nThink you're a wiseguy, huh?  Think you can pull one over on me?\nI solve, and I multiply, and that's it bub.\n\n");
		return;
	}
	
	if(A.rank == 0)
		printf("\n");
	
	PS_free(&A);
	free(x);
	free(y);
}

void run_dot(int argc, char *argv[])
{
	int i,N,na,nd;
	int rank, size;
	double dot;
	double *x, *y;
	
	if(argc < 3)
	{
		printf("Ya need to tell me how big the vectors are kid.  Whadda ya expect me to do, make it up?\n");
		return;
	}
	
	N = atoi(argv[2]);
	
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	na = rank*(size/N);
	if(rank == size-1)
		nd = N - na;
	else
		nd = size/N;
	
	x = (double *) malloc(nd * sizeof(double));
	y = (double *) malloc(nd * sizeof(double));
	
	for(i=0; i<nd; i++)
	{
		x[i] = na + i + 1;
		y[i] = 1.0;
	}
	
	dot = PS_dot(x, y, nd, MPI_COMM_WORLD);
	
	if(rank == 0)
	{
		printf("dot = %lg\nTrue  %lg\n", dot, 0.5*N*(N+1.0));
	}
	
	free(x);
	free(y);
}

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	
	if(argc < 2)
		printf("\nCome on, give me something to work with.\n\n");
	
	else if(strcmp(argv[1], sine_name) == 0)
		run_PS(argc, argv);
	else if(strcmp(argv[1], grav1d_name) == 0)
		run_PS(argc, argv);
	else if(strcmp(argv[1], grav2d_name) == 0)
		run_PS(argc, argv);
	else if(strcmp(argv[1], grav3d_name) == 0)
		run_PS(argc, argv);
	else if(strcmp(argv[1], "dot") == 0)
		run_dot(argc, argv);
	
	else
		printf("\nThat's not even a thing.\n\n");
	
	MPI_Finalize();
	
	return 0;
}

