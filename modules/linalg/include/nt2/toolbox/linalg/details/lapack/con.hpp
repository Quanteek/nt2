/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_DETAILS_LAPACK_CON_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_DETAILS_LAPACK_CON_HPP_INCLUDED
#include <nt2/toolbox/linalg/details/utility/f77_wrapper.hpp>
#include <nt2/toolbox/linalg/details/lapack/workspace.hpp>
// gecon,  pocon


namespace nt2
{
  namespace details
  {
    extern "C"
    {
      #define NT2_COMPLEX void
      void NT2_F77NAME(cgecon)(const char* norm, const long int* n, const NT2_COMPLEX* a,
                           const long int* lda, const float* anorm, float* rcond,
                           NT2_COMPLEX* work, float* rwork, long int* info);
      void NT2_F77NAME(dgecon)(const char* norm, const long int* n, const double* a,
                           const long int* lda, const double* anorm, double* rcond,
                           double* work, long int* iwork, long int* info);
      void NT2_F77NAME(sgecon)(const char* norm, const long int* n, const float* a,
                           const long int* lda, const float* anorm, float* rcond,
                           float* work, long int* iwork, long int* info);
      void NT2_F77NAME(zgecon)(const char* norm, const long int* n, const NT2_COMPLEX* a,
                           const long int* lda, const double* anorm, double* rcond,
                           NT2_COMPLEX* work, double* rwork, long int* info);
      #undef NT2_COMPLEX      
    }
     
#define NT2_GECON(NAME, T)                      \
    inline void gecon(const char* norm,         \
                      const long int* n,        \
                      const T* a,               \
                      const long int* lda,      \
                      const T* anorm,           \
                      T* rcond,                 \
                      long int* info,           \
                      workspace<T> & w)         \
    {                                           \
      w.resizeiw(*n);                           \
      w.resizew(4**n);                          \
      NT2_F77NAME( NAME )(norm, n, a, lda, anorm,                       \
                          rcond, w.getw(), w.getiw(), info);            \
    }                                                                   \
    inline void gecon(const char* norm,                             \
                      const long int* n,                            \
                      const T* a,                                   \
                      const long int* lda,                          \
                      const T* anorm,                               \
                      T* rcond,                                     \
                      long int* info)                               \
    {                                                               \
      workspace<T> w;                                               \
      gecon(norm, n, a, lda, anorm, rcond, info, w);                \
    }                                                               \
        
    NT2_GECON(sgecon, float)
    NT2_GECON(dgecon, double)
      
#undef NT2_GECON
#define NT2_GECON(NAME, T, TBASE)               \
    inline void gecon(const char* norm,         \
                      const long int* n,        \
                      const T* a,               \
                      const long int* lda,      \
                      const TBASE* anorm,       \
                      TBASE* rcond,             \
                      long int* info,           \
                      workspace<T> & w)         \
    {                                           \
      w.resizerw(2**n);                         \
      w.resizew(2**n);                          \
      NT2_F77NAME( NAME )(norm, n, a, lda, anorm,                       \
                          rcond, w.getw(), w.getrw(), info);            \
    }                                                                   \
    inline void gecon(const char* norm,                             \
                      const long int* n,                            \
                      const T* a,                                   \
                      const long int* lda,                          \
                      const TBASE* anorm,                           \
                      TBASE* rcond,                                 \
                      long int* info)                               \
    {                                                               \
      workspace<T> w;                                               \
      gecon(norm, n, a, lda, anorm, rcond, info, w);                \
    }                                                               \
        
    NT2_GECON(cgecon, std::complex<float>, float)
    NT2_GECON(zgecon, std::complex<double>, double)
      
#undef NT2_GECON

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
                               const NT2_COMPLEX* a, const long int* lda,
                               const float* anorm, float* rcond, NT2_COMPLEX* work,
                               float* rwork, long int* info);
      void NT2_F77NAME(zpocon)(const char* uplo, const long int* n,
                               const NT2_COMPLEX* a, const long int* lda,
                               const double* anorm, double* rcond, NT2_COMPLEX* work,
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

#define NT2_POCON(NAME, T)                                              \
    inline void pocon(const char* uplo,                                 \
                      const long int* n,                                \
                      const T* a,                                       \
                      const long int* lda,                              \
                      const T* anorm,                                   \
                      T* rcond,                                         \
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
                      const T* anorm,                                   \
                      T* rcond,                                         \
                          long int* info)                               \
        {                                                               \
          nt2::details::workspace<T> w;                                 \
          pocon(uplo, n, a, lda, anorm, rcond, info, w);                \
        }                                                               \

    NT2_POCON(spocon, float)
    NT2_POCON(dpocon, double)

#undef NT2_POCON    
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
      w.resizerw((*n));                                                 \
      w.resizew(2*(*n));                                                \
      NT2_F77NAME( NAME )(uplo, n, a, lda, anorm,                       \
                      rcond, w.getw(), w.getrw(), info);                \
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

    NT2_POCON(cpocon, std::complex<float>, float)
    NT2_POCON(zpocon, std::complex<double>, double)

#undef NT2_POCON
      }
}


#endif

// /////////////////////////////////////////////////////////////////////////////
// End of con.hpp
// /////////////////////////////////////////////////////////////////////////////
