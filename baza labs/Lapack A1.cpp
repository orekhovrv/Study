#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctime>

/* Auxiliary routines prototypes */
extern "C" {
  extern void dgetrf_(int *M, int *N, double * A, int *LDA, int * IPIV, int *INFO );
  extern void sgetrf_(int *M, int *N, float * A, int *LDA, int * IPIV, int *INFO );
  extern void dgetrs_(char *TRANS, int *N, int *NRHS, double * A, int *LDA,int *IPIV, double *b, int *LDB, int *INFO );
	  
}

/* Main program */
int main() {
        /* Locals */
        int n = 2000;
        int info;
        int ipiv[n];
	double *a; 
	a = (double *)malloc(n*n*sizeof(double));
	double b[n];
	int nrhs = 1;	
	char trans = 'N';
	/*Matrix and b*/
	for ( size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			if (i == j) {
				a[i*n + j] = 1.0;
			} else {
				a[i*n + j] = 1.0/abs(i-j);
			}
		}
		b[i] = sin(double(i));
	}
	/*Calculations*/	
	clock_t start_t = clock();	
	dgetrf_(&n, &n, a, &n, ipiv, &info);
	clock_t end_t = clock();
	double time_rf = double(end_t - start_t)/ double(CLOCKS_PER_SEC);

	start_t = clock();	
	dgetrs_(&trans, &n, &nrhs, a, &n, ipiv, b, &n, &info);
	end_t = clock();
	double time_rs = double(end_t - start_t)/ double(CLOCKS_PER_SEC);	

	printf ("for size %d  init time = %.6f, solve time = %.6f \n ", n, time_rf, time_rs);

       /*for (size_t i=0; i<n; i++) {
        	for (size_t j=0; j<n ; j++)
          printf ("%f ", a[i*n + j]);
        printf("\n");
      }*/
      
      return 0;
} 

