void LU_decomp(int const n, int const lda, double* const A) {

  int i, j, k;

  // Semaphores
  char row_ready[n];
  for (i = 0; i < n; i++)
    row_ready[i] = 0;
  row_ready[0] = 1;

  // For all "iron" rows
  #pragma omp parallel for private(j,k) schedule(dynamic,1)
  for (i = 1; i < n; i++) {
    double * const Ai = A + i*lda; // Pointer to row i

    // For all "hammer" rows
    for (k = 0; k < i; k++) {
      double * const Ak = A + k*lda; // Pointer to row k

      // Spin until "hammer" k is ready
      while (! row_ready[k]) {
        #pragma omp flush
      }

      // Compute the scaling factor (and the element of L)
      Ai[k] /= Ak[k];

      // Hit row "iron" row i with "hammer" row k
      #pragma omp simd
      for (j = k + 1; j < n; j++) 
	Ai[j] -= Ai[k]*Ak[j];

    }

    // Change semaphore for row i to "green"
    row_ready[i] = 1;

  }
}
