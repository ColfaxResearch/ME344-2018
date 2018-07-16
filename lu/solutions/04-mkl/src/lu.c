#include <mkl.h>

void LU_decomp(int const n, int const lda, double* const A) {

  // LU decomposition without pivoting
  // In-place decomposition of form A=LU
  // L is returned below main diagonal of A
  // U is returned at and above main diagonal

  int ipiv[n];
  LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, A, lda, ipiv);

}
