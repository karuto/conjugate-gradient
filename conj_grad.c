/* File:       conj_grad.c
 * Author:     Vincent Zhang
 *
 * Purpose:    A serial conjugate gradient solver program. Due to time limits,
 *             the MPI parallel version is not included in this source code.
 *
 * Compile:    gcc -g -Wall -lm -o conj_grad conj_grad.c
 * Run:        conj_grad [order] [tolerance] [iterations] 
 *                       [Optional suppress output(n)] < [file]
 *
 * Input:      A file that contains a symmetric, positive definite matrix A,  
 *             and the corresponding right hand side vector B. Preferably, each
 *             line consists of [n] elements and the [n+1] line would be the b.
 * Output:     1. The number of iterations,
 *             2. The time used by the solver (not including I/O),
 *             3. The solution to the linear system (if not suppressed),
 *             4. The norm of the residual calculated by the conjugate gradient 
 *                method, and 
 *             5. The norm of the residual calculated directly from the 
 *                definition of residual.
 *
 * Algorithm:  The matrix A's initially read and parsed into an one-dimensional
 *             array; the right hand side vector b is stored in an array as 
 *             well. After some preparation work of allocating memory and 
 *             assigning variables the program jumps into the main loop, the 
 *             conjugate gradient solver. For the exact mathematical procedure,
 *             please refer to http://www.cs.usfca.edu/~peter/cs625/prog2.pdf
 *             and http://en.wikipedia.org/wiki/Conjugate_gradient_method for a
 *             much better demonstration.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
// #include <mpi.h>     /* For MPI functions, etc */
#include "timer.h"


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
	/* The Cross Product */
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
	/* Essentially does the same job as memcpy. */
	int i;
	for (i = 0; i < size; i++) {
		a[i] = b[i];
	}
}




int main(int argc, char **argv) {
  int        i, j, k, order, max_iterations, suppress_output;
  double     start, finish, elapsed;
  double	 tolerance;
   
  // for (i = 0; i < argc; i++) {
  // 	  printf("##### %d = %s\n", i, argv[i]);
  // }
  
  suppress_output = 0;
  
  if (argc > 3) {
	  /* Parse command line arguments */
	  order = atoi(argv[1]); /* Determines dimension of matrix A */
      // printf("==== Dimension of matrix: %d\n", order);
	  
	  tolerance = atof(argv[2]);
	  
	  max_iterations = atoi(argv[3]); 
	  
	  /* Use string compare to check optional suppress */
	  if (argc == 5 && (strcmp(argv[4], "n") == 0)) {
		  printf("==== %s\n", argv[4]);
		  suppress_output = 1;
	  } else {
		  printf("==== No optional command \n");
	  }

	  /* 1-d array for A */
	  double* matrix = malloc(order * order * sizeof(double)); 
 	  /* Right hand side vector, B */ 
	  double* rhs = malloc(order * sizeof(double));
	  
	  
	  /* Read matrix (A) from file and put into array of doubles */
	  for (i = 0; i < order; i++) { /* number of rows */
	  	  for (j = 0; j < order; j++) { /* number of cols */
	        if (!scanf("%lf", &matrix[(i * order) + j])) {
	      		break;
	        }
	      }
	  }
	  
	  /* Read right hand side (B) */
	  for (i = 0; i < order; i++) {
	  	  if (!scanf("%lf", &rhs[i])) {
	    	  break;
	      }
	  }
	  
	  /*
	  printf("==== PRINTING MATRIX ==== \n");
	  for (i = 0; i < (order * order); i++) {
	  	  printf("%f\n", matrix[i]);
	  }
	  printf("==== PRINTING MATRIX ==== \n");
	  */

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
	  
	  GET_TIME(start);
	  
	  while ((k < max_iterations) && (dotProduct(r, r, order) > tolerance)) {
	  	  // memcpy(r_prev_prev, r_prev, (order * sizeof(double)));
		  assignVector(r_prev_prev, r_prev, order);
		  assignVector(r_prev, r, order);
		  assignVector(p_prev, p, order);
	 	  assignVector(x_prev, x, order);
	  	  k++;
		  if (k == 1) {
			  
			  /* P_1 = R_0 */
			  memcpy(p, r_prev, (order * sizeof(double)));
			  memcpy(p_prev, p, (order * sizeof(double)));
			  // assignVector(p, r_prev, order);
			  // assignVector(p_prev, p, order);
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
	  
      GET_TIME(finish);
      elapsed = finish - start;
	  
	  printf("========= Solver Completed ========= \n");
	  printf("The code to be timed took %lf seconds\n", elapsed);
	  printf("Number of iterations: %d \n", k);
	  
	  if (suppress_output == 0) {
		  printf("Solution to the matrix: \n");
		  for (i = 0; i < (order); i++) {
		  	  printf("%f\n", x[i]);
		  }	  	
	  }
	  
	  printf("The norm of the residual calculated by the conjugate gradient method: \n");
	  double norm = 9999.9999;
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
	  


	  /* Free all objects that allocated memory once at the end */
	  free(matrix);
	  free(rhs);
  
	  free(x);
	  free(s);
	  free(p);
	  free(r);
  
	  free(x_prev);
	  free(p_prev);
	  free(r_prev);
	  free(r_prev_prev);
	  free(holderVector);
  
	  
  } else {
	  printf("Usage: %s [order] [tolerance] [iterations] [Optional suppress output(n)]\n", argv[0]);
  }
  
  
  return 0;
} /* main */






