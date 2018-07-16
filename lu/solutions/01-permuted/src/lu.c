void LU_decomp(int const n, int const lda, double* const A) {

  // LU decomposition without pivoting
  // In-place decomposition of form A=LU
  // L is returned below main diagonal of A
  // U is returned at and above main diagonal

  int i, j, k;

  // For all "hammer" rows
  for (k = 0; k < n; k++) {

    double * const Ak = A + k*lda; // Pointer to row k

    // For all "iron" rows
    for (i = k + 1; i < n; i++) {

      double * const Ai = A + i*lda; // Pointer to row i
      
      // Compute the scaling factor (and the element of L)
      Ai[k] /= Ak[k];

      // Hit row "iron" row i with "hammer" row k
      for (j = k + 1; j < n; j++) 
	Ai[j] -= Ai[k]*Ak[j];

    }
  }
}
