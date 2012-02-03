//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II       
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI       
//                                                                              
//          Distributed under the Boost Software License, Version 1.0.          
//                 See accompanying file LICENSE.txt or copy at                 
//                     http://www.boost.org/LICENSE_1_0.txt                     
//==============================================================================
#ifndef NT2_TOOLBOX_ALGEBRA_BLAS_BLAS2_HPP_INCLUDED
#define NT2_TOOLBOX_ALGEBRA_BLAS_BLAS2_HPP_INCLUDED

extern "C"
{
  // Real, double precision
  void F77NAME(dgemv)(const char* trans, const long int* M, const long int* N, 
                      const double* alpha, const double* A, const long int* lda, 
                      const double* dx, const long int* incx, const double* beta, 
                      double* dy, const long int* incy);
  
  void F77NAME(dgbmv)(const char* trans, const long int* M, const long int* N, 
                      const long int* kl, const long int* ku, const double* alpha, 
                      const double* A, const long int* lda, const double* dx, 
                      const long int* incx, const double* beta, double* dy, 
                      const long int* incy);
  
  void F77NAME(dsymv)(const char* uplo, const long int* N, const double* alpha, 
                      const double* A, const long int* lda, const double* dx, 
                      const long int* incx, const double* beta, double* dy, 
                      const long int* incy);
  
  void F77NAME(dsbmv)(const char* uplo, const long int* N, const long int* k, 
                      const double* alpha, const double* A, const long int* lda, 
                      const double* dx, const long int* incx, const double* beta, 
                      double* dy, const long int* incy);
  
  void F77NAME(dspmv)(const char* uplo, const long int* N, const double* alpha,
                      const double* AP, const double* dx, const long int* incx, 
                      const double* beta, double* dy, const long int* incy);
  
  void F77NAME(dtrmv)(const char* uplo, const char* trans, const char* diag,
                      const long int* N, const double* A, const long int* lda, 
                      const double* dx, const long int* incx);
  
  void F77NAME(dtbmv)(const char* uplo, const char* trans, const char* diag, 
                      const long int* N, const long int* K, const double* A, 
                      const long int* lda, double* dx, const long int* incx);
  
  void F77NAME(dtrsv)(const char* uplo, const char* trans, const char* diag, 
                      const long int* N, const double* A, const long int* lda, 
                      double* dx, const long int* incx);

  void F77NAME(dtbsv)(const char* uplo, const char* trans, const char* diag, 
                      const long int* N, const long int* K, const double* A, 
                      const long int* lda, double* dx, const long int* incx);
  
  void F77NAME(dtpsv)(const char* uplo, const char* trans, const char* diag, 
                      const long int* N, double* Ap, double* dx, 
                      const long int* incx);

  // Real, simple precision
  void F77NAME(sger)(const long int* M, const long int* N, const float* alpha, 
                     const float* dx, const long int* incx, const float* dy, 
                     const long int* incy, float* A, const long int* lda);
  
  void F77NAME(ssyr)(const char* uplo, const long int* N, const float* alpha, 
                     const float* dx, const long int* incx, float* A, 
                     const long int* lda);
  
  void F77NAME(sspr)(const char* uplo, const long int* N, const float* alpha, 
                     const float* dx, const long int* incx, float* AP);
  
  void F77NAME(ssyr2)(const char* uplo, const long int* N, const float* alpha, 
                      const float* dx, const long int* incx, const float* dy, 
                      const long int* incy, float* A, const long int* lda);
  
  void F77NAME(sspr2)(const char* uplo, const long int* N, const float* alpha, 
                      const float* dx, const long int* incx, const float* dy, 
                      const long int* incy, float* AP);
  
  void F77NAME(sgemv)(const char* trans, const long int* M, const long int* N, 
                      const float* alpha, const float* A, const long int* lda, 
                      const float* dx, const long int* incx, const float* beta, 
                      float* dy, const long int* incy);
  
  void F77NAME(sgbmv)(const char* trans, const long int* M, const long int* N, 
                      const long int* kl, const long int* ku, const float* alpha, 
                      const float* A, const long int* lda, const float* dx, 
                      const long int* incx, const float* beta, float* dy, 
                      const long int* incy);
  
  void F77NAME(ssymv)(const char* uplo, const long int* N, const float* alpha, 
                      const float* A, const long int* lda, const float* dx, 
                      const long int* incx, const float* beta, float* dy, 
                      const long int* incy);
  
  void F77NAME(ssbmv)(const char* uplo, const long int* N, const long int* k,
                      const float* alpha, const float* A, const long int* lda, 
                      const float* dx, const long int* incx, const float* beta, 
                      float* dy, const long int* incy);
  
  void F77NAME(sspmv)(const char* uplo, const long int* N, const float* alpha, 
                      const float* AP, const float* dx, const long int* incx, 
                      const float* beta, float* dy, const long int* incy);
  
  void F77NAME(strmv)(const char* uplo, const char* trans, const char* diag, 
                      const long int* N, const float* A, const long int* lda, 
                      const float* dx, const long int* incx);

  void F77NAME(stbmv)(const char* uplo, const char* trans, const char* diag, 
                      const long int* N, const long int* K, const float* A, 
                      const long int* lda, const float* dx, const long int* incx);
  
  void F77NAME(strsv)(const char* uplo, const char* trans, const char* diag, 
                      const long int* N, const float* A, float* dx, 
                      const long int* incx);
  
  void F77NAME(stbsv)(const char* uplo, const char* trans, const char* diag, 
                      const long int* N, const float* A, float* dx, 
                      const long int* incx);

  void F77NAME(stbsv)(const char* uplo, const char* trans, const char* diag, 
                      const long int* N, const float* Ap, float* dx, 
                      const long int* incx);
  
  void F77NAME(sger)(const long int* M, const long int* N, const float* alpha, 
                     const float* dx, const long int* incx, const float* dy, 
                     const long int* incy, float* A, const long int* lda);
  
  void F77NAME(ssyr)(const char* uplo, const long int* N, const float* alpha, 
                     const float* dx, const long int* incx, float* A, 
                     const long int* lda);
  
  void F77NAME(sspr)(const char* uplo, const long int* N, const float* alpha, 
                     const float* dx, const long int* incx, float* AP);
  
  void F77NAME(ssyr2)(const char* uplo, const long int* N, const float* alpha, 
                      const float* dx, const long int* incx, const float* dy, 
                      const long int* incy, float* A, const long int* lda);
  
  void F77NAME(sspr2)(const char* uplo, const long int* N, const float* alpha, const float* dx, 
                      const long int* incx, const float* dy, const long int* incy, float* AP);

#endif
