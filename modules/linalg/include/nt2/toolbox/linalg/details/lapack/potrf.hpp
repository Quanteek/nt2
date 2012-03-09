/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_DETAILS_LAPACK_POTRF_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_DETAILS_LAPACK_POTRF_HPP_INCLUDED
#include <nt2/toolbox/linalg/details/utility/f77_wrapper.hpp>
#include <nt2/toolbox/linalg/details/lapack/workspace.hpp>
// The cholesky decomposition is based on potrf

namespace nt2
{
  namespace details
  {

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

  }
}


#endif
