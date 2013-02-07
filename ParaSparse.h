#ifndef PARASPARSE
#define PARASPARSE

#include "mpi.h"

struct ParaSparse
{
	int na;	//First index owned by this rank
	int nd; //Number of indices owned by this rank
	
	int N;	//Global Matrix Size
	int Ne;	//Number of matrix elements
	int *i; //Row Number of entry
	int *j;	//Column Number of entry
	double *Mij; //Entries
	int length;
	
	int rank;	//Process Rank
	int size;	//Number of Processes
	MPI_Comm comm;	//MPI Communicator
	
	int *send_counts;	//Array, size 'size', # of doubles sent to each rank during multiply
	int *recv_counts;	//Array, size 'size', # of doubles received from each rank during multiply
	int *send_displs;	//Array, size 'size', position where rank[i]'s data starts in i_send
	int *recv_displs;	//Array, size 'size', position where rank[i]'s data starts in i_recv
	
	int N_send;		//Number of doubles sent during multiply()
	int N_recv;		//Number of doubles received during multiply()
	int *i_send;	//indices of contributions to y in multiply()
	int *i_recv;	//indices of contributions to y in multiply()
	double *y_send;	//contributions to y in multiply()
	double *y_recv; //contributions to y in multiply()
	int *send_order;
	
	int analyzed;	//flag, whether send_counts is initialized and filled
	int finalized;	//flag, whether i,j,Mij are allocated to the minimuim necessary amount
};

typedef struct ParaSparse ParaSparse;

#define TOL 1e-10

extern const ParaSparse PS_DEFAULT;	

extern char grav1d_name[];
extern char grav2d_name[];
extern char grav3d_name[];
extern char sine_name[];

void PS_generate_sine(ParaSparse *M, int *N, MPI_Comm comm);
void PS_generate_grav1d(ParaSparse *M, int *N, MPI_Comm comm);
void PS_generate_grav2d(ParaSparse *M, int *N, MPI_Comm comm);
void PS_generate_grav3d(ParaSparse *M, int *N, MPI_Comm comm);
void PS_add_entry(ParaSparse *M, int i, int j, double Mij);
void PS_finalize(ParaSparse *M);
void PS_analyze(ParaSparse *M);
void PS_multiply(ParaSparse *M, double *x, double *y);
void PS_printM(ParaSparse *M);
void PS_free(ParaSparse *M);
void PS_printV(double *v, int na, int nd, MPI_Comm comm);
double PS_dot(double *x, double *y, int nd, MPI_Comm comm);
void PS_get_diag(ParaSparse *A, double *d);
void PS_get_precond(ParaSparse *A, ParaSparse *C);
void PS_bcg_iter(ParaSparse *A, ParaSparse *C, double *r, double *p, double *x, 
				 double *rold, double *pold, double *Ap, double *CAp);
int PS_bcg(ParaSparse *A, double *b, double *x);

#endif 
