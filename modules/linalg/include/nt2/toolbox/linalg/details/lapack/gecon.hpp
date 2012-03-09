/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_DETAILS_LAPACK_GECON_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_DETAILS_LAPACK_GECON_HPP_INCLUDED
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
      }
}
#endif

// /////////////////////////////////////////////////////////////////////////////
// End of gecon.hpp
// /////////////////////////////////////////////////////////////////////////////
