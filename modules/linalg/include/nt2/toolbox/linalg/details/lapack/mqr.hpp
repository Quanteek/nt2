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
      void NT2_F77NAME(sormqr)(const char* side, const char* trans,
                               const long int* m, const long int* n, const long int* k,
                               const float* a, const long int* lda,
                               const float* tau, float* c, const long int* ldc,
                               float* work, const long int* lwork, long int* info);
      void NT2_F77NAME(dormqr)(const char* side, const char* trans,
                               const long int* m, const long int* n, const long int* k,
                               const double* a, const long int* lda,
                               const double* tau, double* c, const long int* ldc,
                               double* work, const long int* lwork, long int* info);
      void NT2_F77NAME(zunmqr)(const char* side, const char* trans,
                               const long int* m, const long int* n, const long int* k,
                               const NT2_COMPLEX* a, const long int* lda,
                               const NT2_COMPLEX* tau, NT2_COMPLEX* c, const long int* ldc,
                               NT2_COMPLEX* work, const long int* lwork, long int* info);
      void NT2_F77NAME(cunmqr)(const char* side, const char* trans,
                               const long int* m, const long int* n, const long int* k,
                               const NT2_COMPLEX* a, const long int* lda,
                               const NT2_COMPLEX* tau, NT2_COMPLEX* c, const long int* ldc,
                               NT2_COMPLEX* work, const long int* lwork, long int* info);
#undef NT2_COMPLEX
    }
    
#define NT2_MQR(NAME, T)                        \
   inline void mqr(const char* side,            \
                    const char* trans,          \
                    const long int* m,          \
                    const long int* n,          \
                    const long int* k,          \
                    const T* a,                 \
                    const long int* lda,        \
                    const T* tau,               \
                    T* c,                       \
                    const long int* ldc,        \
                    long int* info,             \
            nt2::details::workspace<T> & w)     \
    {                                           \
      NT2_F77NAME(NAME)(side,trans,             \
                        m,n,k,                  \
                        a,lda,tau,              \
                        c,ldc,                  \
                        w.getw(),w.query(),     \
                        info);                  \
      w.resizew(w.neededsize());                \
      NT2_F77NAME(NAME)(side,trans,             \
                        m,n,k,                  \
                        a,lda,                  \
                        tau,c,ldc,              \
                        w.getw(),&w.neededsize(),\
                        info);                  \
    }                                           \
        inline void mqr(const char* side,       \
                      const char* trans,        \
                      const long int* m,        \
                      const long int* n,        \
                      const long int* k,        \
                      const T* a,               \
                      const long int* lda,      \
                      const T* tau,             \
                      T* c,                     \
                      const long int* ldc,      \
                      long int* info)           \
    {                                           \
      workspace<T> w;                           \
      mqr(side, trans,                          \
            m, n, k,                            \
            a, lda, tau, c, ldc, info, w);      \
    }                                           \

    NT2_MQR(sormqr, float)
    NT2_MQR(dormqr, double)

#undef NT2_MQR

#define NT2_MQR(NAME, T, TBASE)                 \
      inline void mqr(                          \
                      const char* side,         \
                      const char* trans,        \
                      const long int* m,        \
                      const long int* n,        \
                      const long int* k,        \
                      const T* a,               \
                      const long int* lda,      \
                      const T* tau,             \
                      T* c,                     \
                      const long int* ldc,      \
                      long int* info,           \
               nt2::details::workspace<T> & w)  \
      {                                         \
        NT2_F77NAME(NAME)(side,trans,           \
                          m,n,k,a,lda,          \
                          tau,c,ldc,            \
                          w.getw(),w.query(),   \
                          info);                \
        w.resizew(w.neededsize());              \
        NT2_F77NAME(NAME)(side,trans            \
                          ,m,n,k,               \
                          a,lda,                \
                          tau,c,ldc,            \
                      w.getw(),&w.neededsize(), \
                          info);                \
      }                                         \
      inline void mqr(const char* side,         \
                        const char* trans,      \
                        const long int* m,      \
                        const long int* n,      \
                        const long int* k,      \
                        const T* a,             \
                        const long int* lda,    \
                        const T* tau,           \
                        T* c,                   \
                        const long int* ldc,    \
                        long int* info)         \
      {                                         \
        workspace<T> w;                         \
        mqr(side, trans,                        \
              m, n, k,                          \
              a, lda, tau,                      \
              c, ldc, info, w);                 \
      }                                         \

    NT2_MQR(cunmqr, std::complex<float>,  float)
    NT2_MQR(zunmqr, std::complex<double>, double)

#undef NT2_MQR


  }
}


#endif

// /////////////////////////////////////////////////////////////////////////////
// End of qr.hpp
// /////////////////////////////////////////////////////////////////////////////
