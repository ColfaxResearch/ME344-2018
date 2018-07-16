#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>

float IntegrateMyFunction(int const n, float const a, float const b);

void VerifyResult(float const I_numerical, float const a, float const b) {

  // Analytical result
  float const I_exact = 2.0f*(sqrtf(b) - sqrtf(a));

  float const error = fabsf(I_exact - I_numerical)/(fabsf(I_exact) + fabsf(I_numerical));
  if (isnan(I_numerical) || (error > 2e-2f)) {
    printf("ERROR: result of integration is incorrect. I_numerical=%.7e and I_exact=%.7e. Error=%.3e\n", I_numerical, I_exact, error);
    exit(1);
  }

}

double wtime() {
  struct timeval v;
  struct timezone z;

  int i = gettimeofday(&v,&z);
  return ((double)v.tv_sec + (double)v.tv_usec*1.0e-6);
}

int main(const int argc, const char** argv) {

  // Problem size and other parameters
  const int n=10000000;
  const double HztoPerf = 1e-9*double(6*n);

  // Performance measurement
  printf("Numerical integration with n=%d intervals...\n\n", n);

  double rate = 0, dRate = 0; // Benchmarking data
  const int nTrials = 10;
  const int skipTrials = 3; // First few steps are a warm-up
  int trial;
  printf("\033[1m%5s %10s %8s\033[0m\n", "Trial", "Time, s", "GFLOP/s");
  for (trial = 1; trial <= nTrials; trial++) {

    float const a = float(trial);
    float const b = float(trial+1);

    double const tStart = wtime(); // Start timing
    float const I_numerical = IntegrateMyFunction(n, a, b);
    double const tEnd = wtime(); // End timing

    VerifyResult(I_numerical, a, b);

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

}
