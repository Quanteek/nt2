/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2009-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0. 
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_DETAILS_LAPACK_QR_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_DETAILS_LAPACK_QR_HPP_INCLUDED
#include <nt2/toolbox/linalg/details/utility/f77_wrapper.hpp>
#include <nt2/toolbox/linalg/details/lapack/workspace.hpp>
// its a mix or orgqr and ungqr

namespace nt2
{
  namespace details
  {
    extern "C"
    {
      #define NT2_COMPLEX void
      void NT2_F77NAME(dorgqr)(const long int* m, const long int* n, const long int* k,
                           double* a, const long int* lda, const double* tau,
                           double* work, const long int* lwork, long int* info);
      void NT2_F77NAME(sorgqr)(const long int* m, const long int* n, const long int* k,
                           float* a, const long int* lda, const float* tau,
                           float* work, const long int* lwork, long int* info);
      void NT2_F77NAME(zungqr)(const long int* m, const long int* n, const long int* k,
                           NT2_COMPLEX* a, const long int* lda, const NT2_COMPLEX* tau,
                           NT2_COMPLEX* work, const long int* lwork, long int* info);
      void NT2_F77NAME(cungqr)(const long int* m, const long int* n, const long int* k,
                           NT2_COMPLEX* a, const long int* lda, const NT2_COMPLEX* tau,
                           NT2_COMPLEX* work, const long int* lwork, long int* info);
      #undef NT2_COMPLEX
    }
    
#define NT2_GQR(NAME, T)                                                \
    inline void gqr(const long int* m, const long int* n,               \
                    const long int* k,                                  \
                    T* a, const long int* lda,                          \
                    T* tau,                                             \
                    long int* info,                                     \
                    nt2::details::workspace < T >  w)                   \
    {                                                                   \
      NT2_F77NAME( NAME )(m, n, k, a, lda, tau,                             \
                      w.getw(), w.query(), info);                       \
      w.resizew(w.neededsize());                                        \
      NT2_F77NAME( NAME )(m, n, k, a, lda, tau,                             \
                    w.getw(), &w.neededsize(), info);                   \
    }                                                                   \
    inline void gqr(const long int* m, const long int* n,               \
                    const long int* k,                                  \
                    T* a, const long int* lda,                          \
                    T* tau,                                             \
                    long int* info)                                     \
    {                                                                   \
      workspace < T >  w;                                               \
      gqr(m, n, k, a, lda, tau, info, w);                               \
    }                                                                   \
        
    NT2_GQR(dorgqr, double)
    NT2_GQR(sorgqr, float)  
    NT2_GQR(zungqr, std::complex<double>)
    NT2_GQR(cungqr, std::complex<float>)
    
#undef NT2_GQR      
    
    
  }
}


#endif

// /////////////////////////////////////////////////////////////////////////////
// End of qr.hpp
// /////////////////////////////////////////////////////////////////////////////
