/* File:       mpi_float.c
 *
 * Purpose:    A simple Ping-Pong implementation in MPI
 *
 * Compile:    mpicc -g -Wall -o mpi_hello mpi_hello.c
 * Run:        mpiexec -n <number of processes> mpi_hello
 *
 * Input:      None
 * Output:     Length of message sent and average elapsed time per round
 *
 * Algorithm:  It starts with process 0 sending a message to process 1, process
 *             1 receives it and sends 
 *
 * Analysis:   Each process sends a message to process 0,
 *             which prints the messages it has received,
 *             as well as its own message.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
// #include <mpi.h>     /* For MPI functions, etc */
#include "timer.h"

const int MAX_STRING = 100;

int main(int argc, char **argv) {
  char       hostname[MAX_STRING];
  int        my_rank, p, len, i, j, k;
  double     start, finish, elapsed, average_elapsed;
  
  
  double	 tolerance;
  int		 order, max_iterations; 
  char		 suppress[MAX_STRING];
  
  
  // printf("argc = %d\n", argc);
  if (argc == 7 || argc == 8) {
	  
	  // for (int i = 0; i < argc; ++i)
	  // {
	  //     printf("argv[%d]: %s\n", i, argv[i]);
	  // }  	
	  
	  //TODO: Change these as we move into parallel, this is for sequential version only
	  
	  order = atoi(argv[4]); /* Determines dimension of matrix A */
      printf("==== Dimension of matrix: %d\n", order);
	  tolerance = atof(argv[5]);
	  max_iterations = atoi(argv[6]); 
	  if (argc == 8) {
		  printf("==== %s\n", argv[7]);
	  } else {
		  printf("==== No optional command \n");
	  }

	  double* rhs = malloc(order * sizeof(double)); /* Right hand side, b */
	  double* matrix = malloc(order * order * sizeof(double)); /* 1-d array */
	  // for (i = 0; i < order; ++i) { /* number of rows */
	  // 		  matrix[i] = malloc(order*sizeof(double));	/* number of cols */	  	
	  // }
	  FILE *file;
	  file=fopen("matrix.txt", "r");
	  
	  int i, j;
	  
	  /* Read matrix from file and put into array of doubles */
	  for (i = 0; i < order; i++) { /* number of rows */
	  	  for (j = 0; j < order; j++) { /* number of cols */
	      //Use lf format specifier, %c is for character
	        if (!fscanf(file, "%lf", &matrix[(i * order) + j])) {
	      		break;
	        }
		    // printf("%lf\n", matrix[(i * order) + j]); 
		  //Use lf format specifier, \n is for new line
	      }
	  }

	  printf("==== PRINTING MATRIX ==== \n");
	  for (i = 0; i < (order * order); i++) {
	  	  printf("%f\n", matrix[i]);
	  }
	  printf("==== PRINTING MATRIX ==== \n");
	  
	  /* Read right hand side */
	  for (i = 0; i < order; i++) {
	  	  if (!fscanf(file, "%lf", &rhs[i * order])) {
	    	  break;
	      }
		  // printf("%lf\n", rhs[i * order]); 
	  }


	  printf("==== PRINTING RIGHT HAND SIDE ==== \n");
	  for (i = 0; i < order; i++) {
	  	  printf("%f\n", rhs[i]);
	  }
	  printf("==== PRINTING RIGHT HAND SIDE ==== \n");
	  
	  fclose(file);
	  
	  // double k = malloc(sizeof(double));
	  // double x0 = malloc(sizeof(double));
	  
  } else {
	  printf("Incorrect usage. \nUsage: mpiexec -n 4 ./conj_grad 100 1.0e-6 50 (n)\n");
  }
  
  // /* Start up MPI */
  // MPI_Init(NULL, NULL);
  // 
  // /* Get the number of processes */
  // MPI_Comm_size(MPI_COMM_WORLD, &p);
  // 
  // /* Get my rank among all the processes */
  // MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  // 
  // MPI_Get_processor_name(hostname, &len);

  // printf("Proc %d > running on %s\n", my_rank, hostname);
  // 
  // // Allocate space memory for matrices
  // double* A = malloc(1000 * 1000 * sizeof(double));
  // double* B = malloc(1000 * 1000 * sizeof(double));
  // double* C = malloc(1000 * 1000 * sizeof(double));
  // 
  // if (my_rank == 0) {
  //   GET_TIME(start);
  //   
  //   for (i = 0; i < 1000; i++) {
  //     for (j = 0; j < 1000; j++) {
  //       C[i*1000 + j] = 0.0;
  //       for (k = 0; k < 1000; k++) {
  //         C[i*1000 + j] += A[i*1000 + k] * B[k*1000 + j];
  //       }
  //     }
  //   }
  //   
  //   GET_TIME(finish);
  //   elapsed = finish - start;
  //   average_elapsed = elapsed / (2*1000000000);
  //   printf("The code to be timed took %e seconds\n", average_elapsed);
  // }
  // 
  // free(A);

  /* Shut down MPI */
  // MPI_Finalize();
  return 0;
} /* main */


