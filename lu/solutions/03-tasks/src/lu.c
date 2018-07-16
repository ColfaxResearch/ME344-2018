void LU_decomp(int const n, int const lda, double* const A) {

  int i, j, k;

  char row_is_ready[n]; // Semaphores
#pragma omp parallel for
  for (i = 0; i < n; i++)
    row_is_ready[i] = 0;
  row_is_ready[0] = 1;

  #pragma omp parallel private(i, j, k)     // Start up all threads
  {
    #pragma omp master                      // Restrict execution to one thread
    { for (i = 1; i < n; i++) {             // For all "iron" rows
        #pragma omp task firstprivate(i)    // Spawn async task for each "iron"
	{ double * const Ai = A + i*lda;    // Pointer to row i
	  for (k = 0; k < i; k++) { // For all "hammer" rows
	    double * const Ak = A + k*lda; // Pointer to row k
	    // Wait for "hammer" number k to become ready
	    while (!row_is_ready[k]) {
              // Yield to other tasks until we can proceed
	      #pragma omp taskyield
	    }

	    // Compute the scaling factor (and the element of L)
	    Ai[k] /= Ak[k];

	    // Hit row "iron" row i with "hammer" row k
	    #pragma omp simd
	    for (j = k + 1; j < n; j++) 
	      Ai[j] -= Ai[k]*Ak[j];

          }

	  // Done processing the "iron" i. Now it can become a "hammer"
	  row_is_ready[i] = 1;
          #pragma omp flush

	}
      }
    }
  }
}
