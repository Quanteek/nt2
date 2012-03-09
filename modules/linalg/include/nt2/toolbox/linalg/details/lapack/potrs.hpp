/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_DETAILS_LAPACK_POTRS_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_DETAILS_LAPACK_POTRS_HPP_INCLUDED
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
      void NT2_F77NAME(cpotrs)(const char* uplo, const long int* n, const long int* nrhs,
                           const NT2_COMPLEX* a, const long int* lda,
                           NT2_COMPLEX* b, const long int* ldb, long int* info);
      void NT2_F77NAME(zpotrs)(const char* uplo, const long int* n, const long int* nrhs,
                           const NT2_COMPLEX* a, const long int* lda,
                           NT2_COMPLEX* b, const long int* ldb, long int* info);
      void NT2_F77NAME(dpotrs)(const char* uplo, const long int* n, const long int* nrhs,
                           const double* a, const long int* lda,
                           double* b, const long int* ldb, long int* info);
      void NT2_F77NAME(spotrs)(const char* uplo, const long int* n, const long int* nrhs,
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
      NT2_F77NAME( NAME )(uplo, n, nrhs, a, lda,    \
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
