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

double dotProduct(double* a, double* b, int size) {
	double sum = 0;
	int i;
	for (i = 0; i < size; i++) {
		sum += (a[i] * b[i]);
		// printf("Dot product check, round: %d, sum = %lf\n", i, sum);
	}
	// printf("==== Dot product check, final: %lf\n", sum);
	return sum;
}

double* scalarVector(double* dest, double* v, double s, int size) {
	// printf("== Begin: scalar vector product ==");
	int i;
	for (i = 0; i < size; i++) {
		dest[i] = s * v[i];
		// printf("Scalar vector product check, round: %d, sum = %lf\n", i, dest[i]);
	}
	return dest;
}

double* vectorAdd(double* dest, double* a, double* b, int size) {
	// printf("== Begin: vector adding vector ==");
	int i;
	for (i = 0; i < size; i++) {
		dest[i] = a[i] + b[i];
		// printf("Vector add vector check, round: %d, sum = %lf\n", i, dest[i]);
	}
	return dest;
}

double* vectorSubtract(double* dest, double* a, double* b, int size) {
	// printf("== Begin: vector subtracting vector ==");
	int i;
	for (i = 0; i < size; i++) {
		dest[i] = a[i] - b[i];
		// printf("Vector subtract vector check, round: %d, sum = %lf\n", i, dest[i]);
	}
	return dest;
}

double* matrixVector(double* dest, double* matrix, double* v, int size) {
	// printf("== Begin: matrix vector product ==");
	int i, j;
	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			dest[i] += matrix[i*size+j] * v[j];
			// printf("Matrix vector check, round: %d, sum = %lf\n", i, dest[i]);
		}
	}
	return dest;
}





int main(int argc, char **argv) {
  char       hostname[MAX_STRING];
  int        my_rank, p, len, i, j, k;
  double     start, finish, elapsed, average_elapsed;
  
  double	 tolerance;
  int		 order, max_iterations; 
  char		 suppress[MAX_STRING];

  // for (int i = 0; i < argc; ++i)
  // {
  //     printf("argv[%d]: %s\n", i, argv[i]);
  // }  	  
  
  // printf("argc = %d\n", argc);
  if (argc == 7 || argc == 8) {
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
	  
	  int i, j;
	  
	  /* Read matrix from file and put into array of doubles */
	  for (i = 0; i < order; i++) { /* number of rows */
	  	  for (j = 0; j < order; j++) { /* number of cols */
	      //Use lf format specifier, %c is for character
	        if (!scanf("%lf", &matrix[(i * order) + j])) {
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
	  	  if (!scanf("%lf", &rhs[i])) {
	    	  break;
	      }
		  // printf("%lf\n", rhs[i]); 
	  }

	  printf("==== PRINTING RIGHT HAND SIDE ==== \n");
	  for (i = 0; i < order; i++) {
	  	  printf("%f\n", rhs[i]);
	  }
	  printf("==== PRINTING RIGHT HAND SIDE ==== \n");
	  
	  
	  int k = 0;
	  double beta = 0.0;
	  double alpha;
	  double* s = malloc(order * sizeof(double));
	  double* x = malloc(order * sizeof(double));
	  double* x_prev = malloc(order * sizeof(double));
	  double* p = malloc(order * sizeof(double));
	  double* p_prev = malloc(order * sizeof(double));
	  double* r = malloc(order * sizeof(double));
	  double* r_prev = malloc(order * sizeof(double));
	  double* r_prev_prev = malloc(order * sizeof(double));
	  double* tmpVector = malloc(order * sizeof(double));
	  
	  memcpy(r, rhs, (order * sizeof(double)));

	  printf("==== PRINTING R ==== \n");
	  for (i = 0; i < order; i++) {
	  	  printf("%f\n", r[i]);
	  }
	  printf("==== PRINTING R ==== \n");
	  
	  while ((k < max_iterations) && (dotProduct(r, r, order) > tolerance)) {
	      memcpy(x_prev, x, (order * sizeof(double)));
	  	  // memcpy(r_prev_prev, r_prev, (order * sizeof(double)));
	  	  k++;
		  if (k == 1) {
			  r[0] = -8;
			  r[1] = -3;
			  /* P_1 = R_0 */
			  memcpy(p, r, (order * sizeof(double)));
		  } else {
			  /* BETA_k = [R_(k-1) * R_(k-1)] / [R_(k-2) * R_(k-2)]  */
			  beta = dotProduct(r, r, order)/dotProduct(r_prev, r_prev, order);
			  //printf("BETA === %lf\n", beta);-
			  
	      memcpy(p_prev, p, (order * sizeof(double)));
	      
			  /* P_k = R_(k-1) + [BETA_k * P_(k-1)] */
			  tmpVector = scalarVector(tmpVector, p_prev, beta, order);
			  p = vectorAdd(p, r_prev, tmpVector, order);
		  }
		  /* S_k = (A * P_k) */
		  s = matrixVector(s, matrix, p, order);
		  
	  	memcpy(r_prev, r, (order * sizeof(double)));
		  
		  /* ALPHA_k = [R_(k-1) * R_(k-1)] / [P_k * S_k] */
		  double d1 = dotProduct(r_prev, r_prev, order);
		  double d2 = dotProduct(p, s, order);
		  alpha = d1/d2;
		  
		  printf("################## %lf, %lf, %lf\n", d1, d2, alpha);
		  
		  /* X_k = X_(k-1) + (ALPHA_k * P_k) */
		  tmpVector = scalarVector(tmpVector, p, alpha, order);
		  x = vectorAdd(x, x_prev, tmpVector, order);

		  
		  /* R_k = R_(k-1) - (ALPHA_k * S_k) */
		  tmpVector = scalarVector(tmpVector, s, alpha, order);
		  r = vectorSubtract(r, r_prev, tmpVector, order);
	  }
	  printf("========= LOOP COMPLETED =========\n");
	  printf("Number of iterations: %d\n", k);
	  printf("Solution to the matrix:\n");
	  for (i = 0; i < (order); i++) {
	  	  printf("%f\n", x[i]);
	  }
	  printf("The norm of the residual calculated by the conjugate gradient method: \n");
	  for (i = 0; i < (order); i++) {
	  	  printf("%f\n", r[i]);
	  }
	  /*
	  printf("The norm of the residual calculated directly from the definition of residual: \n");
	  //TODO: This is just a holder.
	  for (i = 0; i < (order * order); i++) {
	  	  printf("%f\n", r[i]);
	  }
	  */
	  
	  
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






