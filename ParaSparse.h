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
	MPI_Comm comm;
	
	int *send_counts;
	int *recv_counts;
	int *send_displs;
	int *recv_displs;
	
	int N_send;
	int N_recv;
	int *i_send;
	int *i_recv;
	double *y_send;
	double *y_recv;
	
	int analyzed;
	int finalized;
};

typedef struct ParaSparse ParaSparse;

extern const ParaSparse PS_DEFAULT;	

#define TOL 1e-10

extern char grav1d_name[];
extern char grav2d_name[];
extern char sine_name[];

void PS_generate_sine(ParaSparse *M, int *N, MPI_Comm comm);
void PS_generate_grav1d(ParaSparse *M, int *N, MPI_Comm comm);
void PS_generate_grav2d(ParaSparse *M, int *N, MPI_Comm comm);
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
void PS_bcg_iter(ParaSparse *A, ParaSparse *C, double *r, double *p, double *x);
void PS_bcg(ParaSparse *A, double *b, double *x);

#endif 
