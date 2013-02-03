#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "PS_runner.h"
#include "ParaSparse.h"

void run_sine(int argc, char *argv[])
{
	if(argc < 4)
	{
		printf("Please give problem type (solve/multiply) and a matrix size.\n");
		return;
	}
	
	int i;
	int N = strtod(argv[3], NULL);
	double *x, *y;
	clock_t ticks;
	
	if(N <= 0)
	{
		printf("The size must be positive numbnuts!\n");
		return;
	}
	
	ParaSparse A;
	PS_generate_sine(&A, N, MPI_COMM_WORLD);
	
	x = (double *) malloc(A.nd * sizeof(double));
	y = (double *) malloc(A.nd * sizeof(double));
	
	if(strcmp(argv[2], "multiply") == 0)
	{
		for(i=0; i<A.nd; i++)
			x[i] = A.nd + i+1;
		
		if(A.rank == 0)
			printf("Multiplying a %dx%d system with %d processes.\n\n", N, N, A.size);
		
		PS_analyze(&A);
		MPI_Barrier(MPI_COMM_WORLD);
		if(A.rank == 0)
			ticks = clock();
		
		PS_multiply(&A, x, y);
		
		MPI_Barrier(MPI_COMM_WORLD);
		if(A.rank == 0)
		{
			ticks = clock() - ticks;
			FILE* timefile = fopen("time_sine_mult.out", "a");
			fprintf(timefile, "%d %d %lg\n", N ,(int) ticks, ((double)ticks)/CLOCKS_PER_SEC);
			fclose(timefile);
		}
	}
	
	else if(strcmp(argv[2], "solve") == 0)
	{
		if(A.rank == 0)
			printf("Solving a %dx%d system with %d processes.\n\n", N, N, A.size);
	
		MPI_Barrier(MPI_COMM_WORLD);
		if(A.rank == 0)
			ticks = clock();
		
		PS_bcg(&A, x, y);
		
		MPI_Barrier(MPI_COMM_WORLD);
		if(A.rank == 0)
		{
			ticks = clock() - ticks;
			FILE* timefile = fopen("time_sine_solve.out", "a");
			fprintf(timefile, "%d %d %lg\n", N ,(int) ticks, ((double)ticks)/CLOCKS_PER_SEC);
			fclose(timefile);
		}
	}
	
	PS_free(&A);
	free(x);
	free(y);
}

void run_grav1d(int argc, char *argv[])
{
	if(argc < 4)
	{
		printf("Please give problem type (solve/multiply) and a matrix size.\n");
		return;
	}
	
	int i;
	int N = strtod(argv[3], NULL);
	double *x, *y;
	clock_t ticks;
	
	if(N <= 0)
	{
		printf("The size must be positive numbnuts!\n");
		return;
	}
	
	ParaSparse A;
	PS_generate_grav1d(&A, N, MPI_COMM_WORLD);
	
	x = (double *) malloc(A.nd * sizeof(double));
	y = (double *) malloc(A.nd * sizeof(double));
	
	if(strcmp(argv[2], "multiply") == 0)
	{
		for(i=0; i<A.nd; i++)
			x[i] = A.nd + i+1;
		
		if(A.rank == 0)
			printf("Multiplying a %dx%d system with %d processes.\n\n", N, N, A.size);
		
		PS_analyze(&A);
		MPI_Barrier(MPI_COMM_WORLD);
		if(A.rank == 0)
			ticks = clock();
		
		PS_multiply(&A, x, y);
		
		MPI_Barrier(MPI_COMM_WORLD);
		if(A.rank == 0)
		{
			ticks = clock() - ticks;
			FILE* timefile = fopen("time_grav1d_mult.out", "a");
			fprintf(timefile, "%d %d %lg\n", N ,(int) ticks, ((double)ticks)/CLOCKS_PER_SEC);
			fclose(timefile);
		}
	}
	
	else if(strcmp(argv[2], "solve") == 0)
	{
		if(A.rank == 0)
			printf("Solving a %dx%d system with %d processes.\n\n", N, N, A.size);
		
		MPI_Barrier(MPI_COMM_WORLD);
		if(A.rank == 0)
			ticks = clock();
		
		PS_bcg(&A, x, y);
		
		MPI_Barrier(MPI_COMM_WORLD);
		if(A.rank == 0)
		{
			ticks = clock() - ticks;
			FILE* timefile = fopen("time_grav1d_solve.out", "a");
			fprintf(timefile, "%d %d %lg\n", N ,(int) ticks, ((double)ticks)/CLOCKS_PER_SEC);
			fclose(timefile);
		}
	}
	
	PS_free(&A);
	free(x);
	free(y);
}

void run_grav2d(int argc, char *argv[])
{
	if(argc < 4)
	{
		printf("Please give problem type (solve/multiply) and a number of cells (Nx).\n");
		return;
	}
	
	int i, N;
	int Nx = strtod(argv[3], NULL);
	double *x, *y;
	clock_t ticks;
	
	if(Nx <= 0)
	{
		printf("The size must be positive numbnuts!\n");
		return;
	}
	
	ParaSparse A;
	PS_generate_grav2d(&A, Nx, &N, MPI_COMM_WORLD);
	
	x = (double *) malloc(A.nd * sizeof(double));
	y = (double *) malloc(A.nd * sizeof(double));
	
	if(strcmp(argv[2], "multiply") == 0)
	{
		for(i=0; i<A.nd; i++)
			x[i] = A.nd + i+1;
		
		if(A.rank == 0)
			printf("Multiplying a %dx%d system with %d processes.\n\n", N, N, A.size);
		
		PS_analyze(&A);
		MPI_Barrier(MPI_COMM_WORLD);
		if(A.rank == 0)
			ticks = clock();
		
		PS_multiply(&A, x, y);
		
		MPI_Barrier(MPI_COMM_WORLD);
		if(A.rank == 0)
		{
			ticks = clock() - ticks;
			FILE* timefile = fopen("time_grav2d_mult.out", "a");
			fprintf(timefile, "%d %d %lg\n", N ,(int) ticks, ((double)ticks)/CLOCKS_PER_SEC);
			fclose(timefile);
		}
	}
	
	else if(strcmp(argv[2], "solve") == 0)
	{
		if(A.rank == 0)
			printf("Solving a %dx%d system with %d processes.\n\n", N, N, A.size);
		
		MPI_Barrier(MPI_COMM_WORLD);
		if(A.rank == 0)
			ticks = clock();
		
		PS_bcg(&A, x, y);
		
		MPI_Barrier(MPI_COMM_WORLD);
		if(A.rank == 0)
		{
			ticks = clock() - ticks;
			FILE* timefile = fopen("time_grav2d_solve.out", "a");
			fprintf(timefile, "%d %d %lg\n", N ,(int) ticks, ((double)ticks)/CLOCKS_PER_SEC);
			fclose(timefile);
		}
	}
	
	PS_free(&A);
	free(x);
	free(y);
}

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	
	if(argc < 2)
	{
		printf("Please choose a problem.\n");
		return 0;
	}
	
	else
	{
		if(strcmp(argv[1], sine_name) == 0)
			run_sine(argc, argv);
		
		else if(strcmp(argv[1], grav1d_name) == 0)
			run_grav1d(argc, argv);
		
		else if(strcmp(argv[1], grav2d_name) == 0)
			run_grav2d(argc, argv);
	}
	
	MPI_Finalize();
	
	return 0;
}

