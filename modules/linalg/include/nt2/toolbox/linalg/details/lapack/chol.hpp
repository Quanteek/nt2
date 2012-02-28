/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_DETAILS_LAPACK_CHOL_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_DETAILS_LAPACK_CHOL_HPP_INCLUDED
#include <nt2/toolbox/linalg/details/utility/f77_wrapper.hpp>
#include <nt2/toolbox/linalg/details/lapack/workspace.hpp>
// The cholesky decomposition is based on potrf
// Once done one can call pocon to obtain condition number
// and potrs to solve systems

namespace nt2
{
  namespace details
  {
    //////////////////////////////////////////////////////////////////////
    // pocon calls
    //////////////////////////////////////////////////////////////////////
    /**  purpose         
     **
     **  xpocon estimates the reciprocal of the condition number (in the
     **  1-norm) of a DATA TYPE hermitian positive definite matrix using the
     **  cholesky factorization a = u**h*u or a = l*l**h computed by cpotrf.
     **
     **  an estimate is obtained for norm(inv(a)), and the reciprocal of the
     **  condition number is computed as rcond = 1 / (anorm * norm(inv(a))).
     **
     **  arguments
     **  =========
     **
     **  uplo    (input) char
     **          = 'u':  upper triangle of a is stored;
     **          = 'l':  lower triangle of a is stored.
     **
     **  n       (input) long int
     **          the order of the matrix a.  n >= 0.
     **
     **  a       (input) DATA TYPE array, dimension (lda,n)
     **          the triangular factor u or l from the cholesky factorization
     **          a = u**h*u or a = l*l**h, as computed by cpotrf.
     **
     **  lda     (input) long int
     **          the leading dimension of the array a.  lda >= max(1,n).
     **
     **  anorm   (input) BASE DATA TYPE
     **          the 1-norm (or infinity-norm) of the hermitian matrix a.
     **
     **  rcond   (output) BASE DATA TYPE
     **          the reciprocal of the condition number of the matrix a,
     **          computed as rcond = 1/(anorm * ainvnm), where ainvnm is an
     **          estimate of the 1-norm of inv(a) computed in this routine.
     **
     **
     **
     **  info    (output) long int
     **          = 0:  successful exit
     **          < 0:  if info = -i, the i-th argument had an illegal value
     **/
    extern "C"
    {
      #define NT2_COMPLEX void
      void NT2_F77NAME(cpocon)(const char* uplo, const long int* n,
                               const COMPLEX* a, const long int* lda,
                               const float* anorm, float* rcond, COMPLEX* work,
                               float* rwork, long int* info);
      void NT2_F77NAME(zpocon)(const char* uplo, const long int* n,
                               const COMPLEX* a, const long int* lda,
                               const double* anorm, double* rcond, COMPLEX* work,
                               double* rwork, long int* info);
      void NT2_F77NAME(dpocon)(const char* uplo, const long int* n,
                               const double* a, const long int* lda,
                               const double* anorm, double* rcond, double* work,
                               long int* iwork, long int* info);
      void NT2_F77NAME(spocon)(const char* uplo, const long int* n,
                               const float* a, const long int* lda,
                               const float* anorm, float* rcond, float* work,
                               long int* iwork, long int* info);
      #undef NT2_COMPLEX
    }

#define NT2_POCON(NAME, T, TBASE)                                       \
    inline void pocon(const char* uplo,                                 \
                      const long int* n,                                \
                      const T* a,                                       \
                      const long int* lda,                              \
                      const TBASE* anorm,                               \
                      TBASE* rcond,                                     \
                      long int* info,                                   \
                      workspace<T> & w)                                 \
    {                                                                   \
      w.resizeiw((*n));                                                 \
      w.resizew(3*(*n));                                                \
      NT2_F77NAME( NAME )(uplo, n, a, lda, anorm,                       \
                      rcond, w.getw(), w.getiw(), info);                \
    }                                                                   \
    inline void pocon(const char* uplo,                                 \
                          const long int* n,                            \
                          const T* a,                                   \
                          const long int* lda,                          \
                          const TBASE* anorm,                           \
                          TBASE* rcond,                                 \
                          long int* info)                               \
        {                                                               \
          nt2::details::workspace<T> w;                                 \
          pocon(uplo, n, a, lda, anorm, rcond, info, w);                \
        }                                                               \

    NT2_POCON(spocon, float, float)
    NT2_POCON(dpocon, double, double)
    NT2_POCON(cpocon, std::complex<float>, float)
    NT2_POCON(zpocon, std::complex<double>, double)

#undef NT2_POCON
    //////////////////////////////////////////////////////////////////////
    // potrf calls
    //////////////////////////////////////////////////////////////////////
    /**  purpose
     **  =======
     **
     **  xpotrf computes the cholesky factorization of a DATA TYPE hermitian
     **  positive definite matrix a.
     **
     **  the factorization has the form
     **     a = u**h * u,  if uplo = 'u', or
     **     a = l  * l**h,  if uplo = 'l',
     **  where u is an upper triangular matrix and l is lower triangular.
     **
     **  this is the block version of the algorithm, calling level 3 blas.
     **
     **  arguments
     **  =========
     **
     **  uplo    (input) char
     **          = 'u':  upper triangle of a is stored;
     **          = 'l':  lower triangle of a is stored.
     **
     **  n       (input) long int
     **          the order of the matrix a.  n >= 0.
     **
     **  a       (input/output) DATA TYPE array, dimension (lda,n)
     **          on entry, the hermitian matrix a.  if uplo = 'u', the leading
     **          n-by-n upper triangular part of a contains the upper
     **          triangular part of the matrix a, and the strictly lower
     **          triangular part of a is not referenced.  if uplo = 'l', the
     **          leading n-by-n lower triangular part of a contains the lower
     **          triangular part of the matrix a, and the strictly upper
     **          triangular part of a is not referenced.
     **
     **          on exit, if info = 0, the factor u or l from the cholesky
     **          factorization a = u**h*u or a = l*l**h.
     **
     **  lda     (input) long int
     **          the leading dimension of the array a.  lda >= max(1,n).
     **
     **  info    (output) long int
     **          = 0:  successful exit
     **          < 0:  if info = -i, the i-th argument had an illegal value
     **          > 0:  if info = i, the leading minor of order i is not
     **                positive definite, and the factorization could not be
     **                completed.
     **
     **/      
    extern "C"
    {
      #define NT2_COMPLEX void
      void NT2_F77NAME(dpotrf)(const char* uplo, const long int* n, double* a,
                               const long int* lda, long int* info);
      void NT2_F77NAME(spotrf)(const char* uplo, const long int* n, float* a,
                               const long int* lda, long int* info);
      void NT2_F77NAME(cpotrf)(const char* uplo, const long int* n, NT2_COMPLEX* a,
                               const long int* lda, long int* info);
      void NT2_F77NAME(zpotrf)(const char* uplo, const long int* n, NT2_COMPLEX* a,
                               const long int* lda, long int* info);
      #undef NT2_COMPLEX
    }
    
#define NT2_POTRF(NAME, T)                        \
    inline void potrf(const char* uplo,           \
                      const long int* n,          \
                      T* a,                       \
                      const long int* lda,        \
                      long int* info)             \
    {                                             \
      NT2_F77NAME( NAME )(uplo, n, a, lda, info); \
    }                                             \

    NT2_POTRF(spotrf, float)
    NT2_POTRF(dpotrf, double)
    NT2_POTRF(cpotrf, std::complex<float>)
    NT2_POTRF(zpotrf, std::complex<double>)

#undef NT2_POTRF
      
    //////////////////////////////////////////////////////////////////////
    // potrs calls
    //////////////////////////////////////////////////////////////////////
    /**  purpose
     **  =======
     **
     **  xpotrs solves a system of linear equations a*x = b with a hermitian/symetric
     **  positive definite matrix a using the cholesky factorization 
     **  a = u**(h/t)*u or a = l*l**(h/t) computed by xpotrf.
     **
     **  arguments
     **  =========
     **
     **  uplo    (input) char
     **          = 'u':  upper triangle of a is stored;
     **          = 'l':  lower triangle of a is stored.
     **
     **  n       (input) long int
     **          the order of the matrix a.  n >= 0.
     **
     **  nrhs    (input) long int
     **          the number of right hand sides, i.e., the number of columns
     **          of the matrix b.  nrhs >= 0.
     **
     **  a       (input) DATA TYPE array, dimension (lda,n)
     **          the triangular factor u or l from the cholesky factorization
     **          a = u**h*u or a = l*l**h, as computed by cpotrf.
     **
     **  lda     (input) long int
     **          the leading dimension of the array a.  lda >= max(1,n).
     **
     **  b       (input/output) DATA TYPE array, dimension (ldb,nrhs)
     **          on entry, the right hand side matrix b.
     **          on exit, the solution matrix x.
     **
     **  ldb     (input) long int
     **          the leading dimension of the array b.  ldb >= max(1,n).
     **
     **  info    (output) long int
     **          = 0:  successful exit
     **          < 0:  if info = -i, the i-th argument had an illegal value
     **
     **/
    extern "C"
    {
      #define NT2_COMPLEX void
      void F77NAME(cpotrs)(const char* uplo, const long int* n, const long int* nrhs,
                           const NT2_COMPLEX* a, const long int* lda,
                           NT2_COMPLEX* b, const long int* ldb, long int* info);
      void F77NAME(zpotrs)(const char* uplo, const long int* n, const long int* nrhs,
                           const NT2_COMPLEX* a, const long int* lda,
                           NT2_COMPLEX* b, const long int* ldb, long int* info);
      void F77NAME(dpotrs)(const char* uplo, const long int* n, const long int* nrhs,
                           const double* a, const long int* lda,
                           double* b, const long int* ldb, long int* info);
      void F77NAME(spotrs)(const char* uplo, const long int* n, const long int* nrhs,
                           const float* a, const long int* lda,
                           float* b, const long int* ldb, long int* info);
       #undef NT2_COMPLEX     
    }

#define NT2_POTRS(NAME, T)                      \
    inline void potrs(const char* uplo,         \
                      const long int* n,        \
                      const long int* nrhs,     \
                      const T* a,               \
                      const long int* lda,      \
                      T* b,                     \
                      const long int* ldb,      \
                      long int* info)           \
    {                                           \
      F77NAME( NAME )(uplo, n, nrhs, a, lda,    \
                      b, ldb, info);            \
    }                                           \

    NT2_POTRS(spotrs, float)
    NT2_POTRS(dpotrs, double)
    NT2_POTRS(cpotrs, std::complex<float>)
    NT2_POTRS(zpotrs, std::complex<double>)

#undef NT2_POTRS
    
  }
}


#endif

// /////////////////////////////////////////////////////////////////////////////
// End of l_lse.hpp<3>
// /////////////////////////////////////////////////////////////////////////////
