/* File:       conj_grad.c
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
	double sum = 0.0;
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
		dest[i] = 0.0;
		for (j = 0; j < size; j++) {
			dest[i] += matrix[i*size+j] * v[j];
			// printf("Matrix vector check, round: %d, sum = %lf\n", i, dest[i]);
		}
	}
	return dest;
}

void assignVector(double* a, double* b, int size) {
	int i;
	for (i = 0; i < size; i++) {
		a[i] = b[i];
	}
}



int main(int argc, char **argv) {
  int        p, i, j, k, order, max_iterations;
  double	 tolerance;
   
  // for (i = 0; i < argc; i++) {
  // 	  printf("##### %d = %s\n", i, argv[i]);
  // }
  
  // printf("argc = %d\n", argc);
  if (argc > 3) {
	  order = atoi(argv[1]); /* Determines dimension of matrix A */
      printf("==== Dimension of matrix: %d\n", order);
	  tolerance = atof(argv[2]);
	  max_iterations = atoi(argv[3]); 
	  if (argc == 5) {
		  printf("==== %s\n", argv[4]);
	  } else {
		  printf("==== No optional command \n");
	  }

	  /* 1-d array for A */
	  double* matrix = malloc(order * order * sizeof(double)); 
 	  /* Right hand side, b */ 
	  double* rhs = malloc(order * sizeof(double));
	  
	  
	  /* Read matrix (A) from file and put into array of doubles */
	  for (i = 0; i < order; i++) { /* number of rows */
	  	  for (j = 0; j < order; j++) { /* number of cols */
	        if (!scanf("%lf", &matrix[(i * order) + j])) {
	      		break;
	        }
	      }
	  }
	  
	  /* Read right hand side (b) */
	  for (i = 0; i < order; i++) {
	  	  if (!scanf("%lf", &rhs[i])) {
	    	  break;
	      }
	  }
	  
	  
	  printf("==== PRINTING MATRIX ==== \n");
	  for (i = 0; i < (order * order); i++) {
	  	  printf("%f\n", matrix[i]);
	  }
	  printf("==== PRINTING MATRIX ==== \n");
	  
	  

	  /*
	  printf("==== PRINTING RIGHT HAND SIDE ==== \n");
	  for (i = 0; i < order; i++) {
	  	  printf("%f\n", rhs[i]);
	  }
	  printf("==== PRINTING RIGHT HAND SIDE ==== \n");
	  */
	  
	  
	  /* Prepare variables for the main loop */
	  k = 0;
	  double beta;
	  double alpha;
	  double* x = malloc(order * sizeof(double));
	  double* s = malloc(order * sizeof(double));
	  double* x_prev = malloc(order * sizeof(double));
	  double* p = malloc(order * sizeof(double));
	  double* p_prev = malloc(order * sizeof(double));
	  double* r = malloc(order * sizeof(double));
	  double* r_prev = malloc(order * sizeof(double));
	  double* r_prev_prev = malloc(order * sizeof(double));
	  double* holderVector = malloc(order * sizeof(double));
	  
	  /* Before we start, copy values of B into R_0 */
	  // memcpy(r, rhs, (order * sizeof(double)));
	  // memcpy(r_prev, rhs, (order * sizeof(double)));
	  // memcpy(r_prev_prev, rhs, (order * sizeof(double)));
	  for (i = 0; i < order; i++) {
		  x[i] = 0;
		  r[i] = rhs[i];
		  r_prev[i] = rhs[i];
		  r_prev_prev[i] = rhs[i];
	  }
	  
	  while ((k < max_iterations) && (dotProduct(r, r, order) > tolerance)) {
	  	  // memcpy(r_prev_prev, r_prev, (order * sizeof(double)));
		  assignVector(r_prev_prev, r_prev, order);
		  assignVector(r_prev, r, order);
		  assignVector(p_prev, p, order);
	 	  assignVector(x_prev, x, order);
	  	  k++;
		  if (k == 1) {
			  
			  /* P_1 = R_0 */
			  // memcpy(p, r_prev, (order * sizeof(double)));
			  // memcpy(p_prev, p, (order * sizeof(double)));
			  assignVector(p, r_prev, order);
			  assignVector(p_prev, p, order);
		  } else {
			  /* BETA_k = [R_(k-1) * R_(k-1)] / [R_(k-2) * R_(k-2)]  */
			  beta = dotProduct(r_prev, r_prev, order)/dotProduct(r_prev_prev, r_prev_prev, order);
			  
	          // memcpy(p_prev, p, (order * sizeof(double)));
	      
			  /* P_k = R_(k-1) + [BETA_k * P_(k-1)] */
			  holderVector = scalarVector(holderVector, p_prev, beta, order);
			  p = vectorAdd(p, r, holderVector, order);
		  }
		  /* S_k = (A * P_k) */
		  s = matrixVector(s, matrix, p, order);
	  
	  	  // memcpy(r_prev, r, (order * sizeof(double)));
		  
		  /* ALPHA_k = [R_(k-1) * R_(k-1)] / [P_k * S_k] */
		  double d1 = dotProduct(r_prev, r_prev, order);
		  double d2 = dotProduct(p, s, order);
		  alpha = d1/d2;
		  
	      // memcpy(x_prev, x, (order * sizeof(double)));
	  		  
		  /* X_k = X_(k-1) + (ALPHA_k * P_k) */
		  holderVector = scalarVector(holderVector, p, alpha, order);
		  x = vectorAdd(x, x_prev, holderVector, order);
		  
		  /* R_k = R_(k-1) - (ALPHA_k * S_k) */
		  holderVector = scalarVector(holderVector, s, alpha, order);
		  r = vectorSubtract(r, r_prev, holderVector, order);
		  
	  }
	  
	  printf("========= Solver Completed =========\n");
	  printf("Number of iterations: %d\n", k);
	  printf("Solution to the matrix:\n");
	  for (i = 0; i < (order); i++) {
	  	  printf("%f\n", x[i]);
	  }
	  
	  printf("The norm of the residual calculated by the conjugate gradient method: \n");
	  double norm;
	  for (i = 0; i < (order); i++) {
	  	  norm += r[i]*r[i];
	  }
	  printf("%lf\n", sqrt(norm));
	  
	  
	  /* Calculate the residual with the algorithm
	  	 R_k = B - (A * X_k) */
	  r = vectorSubtract(r, rhs, matrixVector(holderVector, matrix, x, order), order);
	  
	  
	  printf("The norm of the residual calculated directly from the definition of residual: \n");
	  for (i = 0; i < (order * order); i++) {
	  	  norm += r[i]*r[i];
	  }
	  printf("%lf\n", sqrt(norm));
	  
	  
	  
  } else {
	  printf("Usage: %s [order] [tolerance] [iterations] [Suppress(n)]\n", argv[0]);
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






