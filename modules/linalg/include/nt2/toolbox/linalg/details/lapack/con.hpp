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



namespace nt2
{
  namespace details
  {
    extern "C"
    {
      #define NT2_COMPLEX void
      void F77NAME(cgecon)(const char* norm, const long int* n, const COMPLEX* a,
                           const long int* lda, const float* anorm, float* rcond,
                           COMPLEX* work, float* rwork, long int* info);
      void F77NAME(dgecon)(const char* norm, const long int* n, const double* a,
                           const long int* lda, const double* anorm, double* rcond,
                           double* work, long int* iwork, long int* info);
      void F77NAME(sgecon)(const char* norm, const long int* n, const float* a,
                           const long int* lda, const float* anorm, float* rcond,
                           float* work, long int* iwork, long int* info);
      void F77NAME(zgecon)(const char* norm, const long int* n, const COMPLEX* a,
                           const long int* lda, const double* anorm, double* rcond,
                           COMPLEX* work, double* rwork, long int* info);
      #undef NT2_COMPLEX      
    }
    
#define LPP_GECON(NAME, T, TBASE)               \
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
      F77NAME( NAME )(norm, n, a, lda, anorm,                           \
                      rcond, w.getw(), w.getrw(), info);                \
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
        
    LPP_GECON(sgecon, float, float)
    LPP_GECON(dgecon, double, double)
    LPP_GECON(cgecon, std::complex<float>, float)
    LPP_GECON(zgecon, std::complex<double>, double)
      
#undef LPP_GECON

    
  }
}


#endif

// /////////////////////////////////////////////////////////////////////////////
// End of con.hpp
// /////////////////////////////////////////////////////////////////////////////
