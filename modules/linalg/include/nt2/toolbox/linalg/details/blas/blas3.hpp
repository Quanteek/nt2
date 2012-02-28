//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II       
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI       
//                                                                              
//          Distributed under the Boost Software License, Version 1.0.          
//                 See accompanying file LICENSE.txt or copy at                 
//                     http://www.boost.org/LICENSE_1_0.txt                     
//==============================================================================
#ifndef NT2_TOOLBOX_LINALG_DETAILS_BLAS_BLAS3_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_DETAILS_BLAS_BLAS3_HPP_INCLUDED
#include <nt2/toolbox/linalg/details/utility/f77_wrapper.hpp>

extern "C"
{
// Real, single precision
  void NT2_F77NAME(sgemm)(const char *transa, const char *transb, const long int *m, 
                      const long int *n, const long int *k, const float *alpha, 
                      const float *a, const long int *lda, const float *b, 
                      const long int *ldb, const float *beta, float *c, 
                      const long int *ldc);
  
  void NT2_F77NAME(strsm)(const char *side, const char *uplo, const char *transa, 
                      const char *diag, const long int *m, const long int *n, 
                      float *alpha, const float *A, const long int *lda,
                      const float *B, const long int *ldb);
  
  void NT2_F77NAME(strmm)(const char *side, const char *uplo, const char *transa, 
                      const char *diag, const long int *m, const long int *n, 
                      float *alpha, const float *A, const long int *lda,
                      float *B, const long int *ldb);
  
  void NT2_F77NAME(ssymm)(const char *side, const char *uplo, const long int *m, 
                      const long int *n, const float *alpha, const float *A, 
                      const long int *lda, const float *B, const long int *ldb, 
                      const float *beta, float *C, const long int *ldc);
  
  void NT2_F77NAME(ssyrk)(const char *uplo, const char *transa, const long int *n, 
                      const long int *k, const float *alpha, const float *A, 
                      const long int *lda, const float *beta, float *C, 
                      const long int *ldc);
  
  void NT2_F77NAME(ssyr2k)(const char *uplo, const char *transa, const long int *n, 
                       const long int *k, const float *alpha, const float *A, 
                       const long int *lda, const float *B, const long int *ldb,
                       const float *beta, float *C, const long int *ldc);


// Real, double precison
  void NT2_F77NAME(dgemm)(const char *transa, const char *transb, const long int *m, 
                      const long int *n, const long int *k, const double *alpha, 
                      const double *a, const long int *lda, const double *b, 
                      const long int *ldb, const double *beta, double *c, 
                      const long int *ldc);
  
  void NT2_F77NAME(dtrsm)(const char *side, const char *uplo, const char *transa, 
                      const char *diag, const long int *m, const long int *n, 
                      const double *alpha, const double *A, const long int *lda,
                      double *B, const long int *ldb);
  
  void NT2_F77NAME(dtrmm)(const char *side, const char *uplo, const char *transa, 
                      const char *diag, const long int *m, const long int *n, 
                      const double *alpha, const double *A, const long int *lda,
                      double *B, const long int *ldb);
  
  void NT2_F77NAME(dsymm)(const char *side, const char *uplo, const long int *m, 
                      const long int *n, const double *alpha, const double *A, 
                      const long int *lda, const double *B, 
                      const long int *ldb, const double *beta, double *C, 
                      const long int *ldc);
  
  void NT2_F77NAME(dsyrk)(const char *uplo, const char *transa, const long int *n, 
                      const long int *k, const double *alpha, const double *A, 
                      const long int *lda, const double *beta, double *C, 
                      const long int *ldc);
  
  void NT2_F77NAME(dsyr2k)(const char *uplo, const char *transa, const long int *n, 
                       const long int *k, const double *alpha, const double *A, 
                       const long int *lda, const double *B, const long int *ldb,
                       const double *beta, double *C, const long int *ldc);

// Complex, single precision
#define NT2_WRAP_COMPLEX void

  void NT2_F77NAME(cgemm)(const char *transa, const char *transb, const long int *m, 
                      const long int *n, const long int *k,
                      const NT2_WRAP_COMPLEX *alpha, const NT2_WRAP_COMPLEX *a, 
                      const long int *lda, const NT2_WRAP_COMPLEX *b, 
                      const long int *ldb, const NT2_WRAP_COMPLEX *beta, 
                      NT2_WRAP_COMPLEX *c, const long int *ldc);
  
  void NT2_F77NAME(ctrsm)(const char *side, const char *uplo, const char *transa, 
                      const char *diag, const long int *m, const long int *n, 
                      const NT2_WRAP_COMPLEX *alpha, const NT2_WRAP_COMPLEX *A, 
                      const long int *lda, NT2_WRAP_COMPLEX *B, 
                      const long int *ldb);
  
  void NT2_F77NAME(ctrmm)(const char *side, const char *uplo, const char *transa, 
                      const char *diag, const long int *m, const long int *n, 
                      NT2_WRAP_COMPLEX *alpha, const NT2_WRAP_COMPLEX *A, 
                      const long int *lda, const NT2_WRAP_COMPLEX *B, 
                      const long int *ldb);
  
  void NT2_F77NAME(csymm)(const char *side, const char *uplo, const long int *m, 
                      const long int *n, const NT2_WRAP_COMPLEX *alpha, 
                      const NT2_WRAP_COMPLEX *A, const long int *lda, 
                      const NT2_WRAP_COMPLEX *B, const long int *ldb, 
                      const NT2_WRAP_COMPLEX *beta, NT2_WRAP_COMPLEX *C, 
                      const long int *ldc);

  void NT2_F77NAME(chemm)(const char *side, const char *uplo, const long int *m, 
                      const long int *n, const NT2_WRAP_COMPLEX *alpha, 
                      const NT2_WRAP_COMPLEX *A, const long int *lda, 
                      const NT2_WRAP_COMPLEX *B, const long int *ldb, 
                      const NT2_WRAP_COMPLEX *beta, NT2_WRAP_COMPLEX *C, 
                      const long int *ldc);
  
  void NT2_F77NAME(csyrk)(const char *uplo, const char *transa, const long int *n, 
                      const long int *k, const NT2_WRAP_COMPLEX *alpha, 
                      const NT2_WRAP_COMPLEX *A, const long int *lda, 
                      const NT2_WRAP_COMPLEX *beta, NT2_WRAP_COMPLEX *C, 
                      const long int *ldc);
  
  void NT2_F77NAME(csyr2k)(const char *uplo, const char *transa, const long int *n, 
                       const long int *k, const NT2_WRAP_COMPLEX *alpha, 
                       const NT2_WRAP_COMPLEX *A, const long int *lda, 
                       const NT2_WRAP_COMPLEX *B, const long int *ldb,
                       const NT2_WRAP_COMPLEX *beta, NT2_WRAP_COMPLEX *C, 
                       const long int *ldc);

// Complex, double precision
  void NT2_F77NAME(zgemm)(const char *transa, const char *transb, const long int *m, 
                      const long int *n, const long int *k,
                      const NT2_WRAP_COMPLEX *alpha, const NT2_WRAP_COMPLEX *a, 
                      const long int *lda, const NT2_WRAP_COMPLEX *b, 
                      const long int *ldb, const NT2_WRAP_COMPLEX *beta, 
                      NT2_WRAP_COMPLEX *c, const long int *ldc);
  
  void NT2_F77NAME(ztrsm)(const char *side, const char *uplo, const char *transa, 
                      const char *diag, const long int *m, const long int *n, 
                      const NT2_WRAP_COMPLEX *alpha, const NT2_WRAP_COMPLEX *A, 
                      const long int *lda, const NT2_WRAP_COMPLEX *B, 
                      const long int *ldb);
  
  void NT2_F77NAME(ztrmm)(const char *side, const char *uplo, const char *transa, 
                      const char *diag, const long int *m, const long int *n, 
                      const NT2_WRAP_COMPLEX *alpha, const NT2_WRAP_COMPLEX *A, 
                      const long int *lda, const NT2_WRAP_COMPLEX *B, 
                      const long int *ldb);
  
  void NT2_F77NAME(zhemm)(const char *side, const char *uplo, const long int *m, 
                      const long int *n, const NT2_WRAP_COMPLEX *alpha, 
                      const NT2_WRAP_COMPLEX *A, const long int *lda, 
                      const NT2_WRAP_COMPLEX *B, const long int *ldb, 
                      const NT2_WRAP_COMPLEX *beta, NT2_WRAP_COMPLEX *C, 
                      const long int *ldc);

  void NT2_F77NAME(zsymm)(const char *side, const char *uplo, const long int *m, 
                      const long int *n, const NT2_WRAP_COMPLEX *alpha, 
                      const NT2_WRAP_COMPLEX *A, const long int *lda,
                      const NT2_WRAP_COMPLEX *B, const long int *ldb, 
                      const NT2_WRAP_COMPLEX *beta, NT2_WRAP_COMPLEX *C, 
                      const long int *ldc);
  
  void NT2_F77NAME(zsyrk)(const char *uplo, const char *transa, const long int *n, 
                      const long int *k, const NT2_WRAP_COMPLEX *alpha, 
                      const NT2_WRAP_COMPLEX *A, const long int *lda, 
                      const NT2_WRAP_COMPLEX *beta, NT2_WRAP_COMPLEX *C, 
                      const long int *ldc);
  
  void NT2_F77NAME(zsyr2k)(const char *uplo, const char *transa, const long int *n, 
                       const long int *k, const NT2_WRAP_COMPLEX *alpha, 
                       const NT2_WRAP_COMPLEX *A, const long int *lda, 
                       const NT2_WRAP_COMPLEX *B, const long int *ldb,
                       const NT2_WRAP_COMPLEX *beta, NT2_WRAP_COMPLEX *C, 
                       const long int *ldc);
}

#endif
