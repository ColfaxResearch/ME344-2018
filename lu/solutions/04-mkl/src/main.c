#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>

void LU_decomp(const int n, const int lda, double* const A);

void VerifyResult(const int n, const int lda, double* LU, double* refA) {

  // Verifying that A=LU
  double A[n*lda];
  double L[n*lda];
  double U[n*lda];
  int c, i, j, k;
  for (c = 0; c < n*lda; c++)
    A[c] = L[c] = U[c] = 0.0;

  for (i = 0; i < n; i++) {
    for (j = 0; j < i; j++)
      L[i*lda + j] = LU[i*lda + j];
    L[i*lda+i] = 1.0;
    for (j = i; j < n; j++)
      U[i*lda + j] = LU[i*lda + j];
  }
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      for (k = 0; k < n; k++)
	A[i*lda + j] += L[i*lda + k]*U[k*lda + j];

  double deviation1 = 0.0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      deviation1 += (refA[i*lda+j] - A[i*lda+j])*(refA[i*lda+j] - A[i*lda+j]);
    }
  }
  deviation1 /= (double)(n*lda);
  if (isnan(deviation1) || (deviation1 > 1.0e-2)) {
    printf("ERROR: LU is not equal to A (deviation1=%e)!\n", deviation1);
    exit(1);
  }

#ifdef VERBOSE
  printf("\n(L-D)+U:\n");
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++)
      printf("%10.3e", LU[i*lda+j]);
    printf("\n");
  }

  printf("\nL:\n");
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++)
      printf("%10.3e", L[i*lda+j]);
    printf("\n");
  }

  printf("\nU:\n");
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++)
      printf("%10.3e", U[i*lda+j]);
    printf("\n");
  }

  printf("\nLU:\n");
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++)
      printf("%10.3e", A[i*lda+j]);
    printf("\n");
  }

  printf("\nA:\n");
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++)
      printf("%10.3e", refA[i*lda+j]);
    printf("\n");
  }

  printf("deviation1=%e\n", deviation1);
#endif

}

double wtime() {
  struct timeval v;
  struct timezone z;

  int i = gettimeofday(&v,&z);
  return ((double)v.tv_sec + (double)v.tv_usec*1.0e-6);
}

int main(const int argc, const char** argv) {

  // Problem size and other parameters
  const int n=512;
  const int lda=528;
  const double HztoPerf = 1e-9*2.0/3.0*double(n*n*lda);

  const size_t containerSize = sizeof(double)*n*lda+64;
  double* A = (double*) _mm_malloc(containerSize, 64);
  double* B = (double*) _mm_malloc(containerSize, 64);

  int i, j, c;

  // Initialize matrix
  for (i = 0; i < n; i++) {
    double sum = 0.0;
    for (j = 0; j < n; j++) {
      A[i*lda+j] = (double)(i*n+j);
      sum += A[i*lda+j];
    }
    sum -= A[i*lda+i];
    A[i*lda+i] = 2.0f*sum;
  }
  A[(n-1)*lda+n] = 0.0; // Touch just in case

  for (c = 0; c < n*lda; c++)
    B[c] = A[c];
  
  // Performance measurement
  printf("LU decomposition of a diagonally-dominant matrix of size %dx%d...\n\n", n, n	 );

  double rate = 0, dRate = 0; // Benchmarking data
  const int nTrials = 10;
  const int skipTrials = 3; // First step is warm-up on Xeon Phi coprocessor
  int trial;
  printf("\033[1m%5s %10s %8s\033[0m\n", "Trial", "Time, s", "GFLOP/s");
  for (trial = 1; trial <= nTrials; trial++) {

    for (c = 0; c < n*lda; c++)
      A[c] = B[c];

    const double tStart = wtime(); // Start timing
    LU_decomp(n, lda, A);
    const double tEnd = wtime(); // End timing

    VerifyResult(n, lda, A, B);

    if (trial > skipTrials) { // Collect statistics
      rate  += HztoPerf/(tEnd - tStart); 
      dRate += HztoPerf*HztoPerf/((tEnd - tStart)*(tEnd-tStart)); 
    }

    printf("%5d %10.3e %8.2f %s\n", 
	   trial, (tEnd-tStart), HztoPerf/(tEnd-tStart), (trial<=skipTrials?"*":""));
    fflush(stdout);
  }
  rate/=(double)(nTrials-skipTrials); 
  dRate=sqrt(dRate/(double)(nTrials-skipTrials)-rate*rate);
  printf("-----------------------------------------------------\n");
  printf("\033[1m%s %4s \033[42m%10.2f +- %.2f GFLOP/s\033[0m\n",
	 "Average performance:", "", rate, dRate);
  printf("-----------------------------------------------------\n");
  printf("* - warm-up, not included in average\n\n");

  _mm_free(A);
  _mm_free(B);

}
