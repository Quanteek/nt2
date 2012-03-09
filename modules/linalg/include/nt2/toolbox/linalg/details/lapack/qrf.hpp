/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2009-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI 
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_DETAILS_LAPACK_QRF_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_DETAILS_LAPACK_QRF_HPP_INCLUDED
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
      void NT2_F77NAME(dgeqrf)(const long int* m, const long int* n,
                           double* a, const long int* lda, double* tau,
                           double* work, const long int* lwork, long int* info);
      void NT2_F77NAME(sgeqrf)(const long int* m, const long int* n,
                           float* a, const long int* lda, float* tau,
                           float* work, const long int* lwork, long int* info);
      void NT2_F77NAME(zgeqrf)(const long int* m, const long int* n,
                           NT2_COMPLEX* a, const long int* lda, NT2_COMPLEX* tau,
                           NT2_COMPLEX* work, const long int* lwork, long int* info);
      void NT2_F77NAME(cgeqrf)(const long int* m, const long int* n,
                           NT2_COMPLEX* a, const long int* lda, NT2_COMPLEX* tau,
                           NT2_COMPLEX* work, const long int* lwork, long int* info);
      #undef NT2_COMPLEX
    }
    
#define NT2_GEQRF(NAME, T)                      \
    inline void geqrf(const long int* m,        \
                      const long int* n,        \
                      T* a,                     \
                      const long int* lda,      \
                      T* tau,                   \
                      long int* info,           \
                nt2::details::workspace<T> & w) \
    {                                           \
      NT2_F77NAME( NAME )(m, n, a, lda, tau,    \
                      w.getw(),w.query(),info); \
      w.resizew(w.neededsize());                \
      NT2_F77NAME( NAME )(m, n, a, lda, tau,    \
                      w.getw(),                 \
                      &w.neededsize(),info);    \
    }                                           \
    inline void geqrf(const long int* m,        \
                      const long int* n,        \
                      T* a,                     \
                      const long int* lda,      \
                      T* tau,                   \
                      long int* info)           \
    {                                           \
      nt2::details::workspace<T> w;             \
      geqrf(m, n, a, lda, tau, info, w);        \
    }                                           \

    NT2_GEQRF(sgeqrf, float)
    NT2_GEQRF(dgeqrf, double)
    NT2_GEQRF(cgeqrf, std::complex<float>)
    NT2_GEQRF(zgeqrf, std::complex<double>)

#undef NT2_GEQRF
    
  }
}


#endif

// /////////////////////////////////////////////////////////////////////////////
// End of qrf.hpp
// /////////////////////////////////////////////////////////////////////////////
