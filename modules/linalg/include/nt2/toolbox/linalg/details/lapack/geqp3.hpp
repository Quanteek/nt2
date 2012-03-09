/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P 
 *         Copyright 2009-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_DETAILS_LAPACK_GEQP3_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_DETAILS_LAPACK_GEQP3_HPP_INCLUDED
#include <nt2/toolbox/linalg/details/utility/f77_wrapper.hpp>
#include <nt2/toolbox/linalg/details/lapack/workspace.hpp>
// its geqp3

namespace nt2
{
  namespace details
  {
    extern "C"
    {
#define NT2_COMPLEX void
      void NT2_F77NAME(cgeqp3)(const long int* m, const long int* n,
                               NT2_COMPLEX* a, const long int* lda,
                               long int* jpvt, NT2_COMPLEX* tau, NT2_COMPLEX* work,
                               const long int* lwork, float* rwork, long int* info);
      void NT2_F77NAME(dgeqp3)(const long int* m, const long int* n,
                               double* a, const long int* lda,
                               long int* jpvt, double* tau, double* work,
                               const long int* lwork, long int* info);
      void NT2_F77NAME(sgeqp3)(const long int* m, const long int* n,
                               float* a, const long int* lda,
                               long int* jpvt, float* tau, float* work,
                               const long int* lwork, long int* info);
      void NT2_F77NAME(zgeqp3)(const long int* m, const long int* n,
                               NT2_COMPLEX* a, const long int* lda,
                               long int* jpvt, NT2_COMPLEX* tau, NT2_COMPLEX* work,
                               const long int* lwork, double* rwork, long int* info);
#undef NT2_COMPLEX
    }
    

#define NT2_GEQP3(NAME, T)                              \
    inline void geqp3(const long int* m,                \
                      const long int* n,                \
                      T* a,                             \
                      const long int* lda,              \
                      long int* jpvt,                   \
                      T* tau,                           \
                      long int* info,                   \
                      nt2::details::workspace<T> & w)   \
    {                                                   \
      NT2_F77NAME(NAME)(m,n,a,lda,jpvt,tau,             \
                        w.getw(),w.query(),info);       \
      w.resizew(w.neededsize());                        \
      NT2_F77NAME(NAME)(m,n,a,lda,jpvt,tau,             \
                        w.getw(),&w.neededsize(),info); \
    }                                                   \
    inline void geqp3(const long int* m,                \
                      const long int* n,                \
                      T* a,                             \
                      const long int* lda,              \
                      long int* jpvt,                   \
                      T* tau,                           \
                      long int* info)                   \
    {                                                   \
      nt2::details::workspace<T> w;                     \
      geqp3(m, n, a, lda, jpvt,                         \
            tau, info, w);                              \
    }                                                   \
        
    NT2_GEQP3(sgeqp3, float)
    NT2_GEQP3(dgeqp3, double)

#undef NT2_GEQP3

#define NT2_GEQP3(NAME,T,TBASE)                 \
    inline void geqp3(const long int *m,        \
                      const long int *n,        \
                      T*a,                      \
                      const long int *lda,      \
                      long int* jpvt,           \
                      T*tau,                    \
                      long int *info,           \
              nt2::details::workspace<T>&w)     \
      {                                         \
        w.resizerw(2**n);                       \
        NT2_F77NAME(NAME)(m,n,a,lda,jpvt,tau,   \
                      w.getw(),w.query(),       \
                      w.getrw(),info);          \
        w.resizew(w.neededsize());              \
        NT2_F77NAME(NAME)(m,n,a,lda,jpvt,tau,   \
                w.getw(),&w.neededsize(),       \
                w.getrw(),info);                \
  }                                             \
  inline void geqp3(const long int *m,          \
                  const long int *n,            \
                  T*a,                          \
                  const long int *lda,          \
                  long int * jpvt,              \
                  T* tau,                       \
                  long int *info)               \
  {                                             \
    nt2::details::workspace<T>w;                \
    geqp3(m,n,a,lda,jpvt,tau,info,w);           \
  }                                             \

    NT2_GEQP3(cgeqp3, std::complex<float>,  float)
    NT2_GEQP3(zgeqp3, std::complex<double>, double)

#undef NT2_GEQP3

    
  }
}


#endif

// /////////////////////////////////////////////////////////////////////////////
// End of geqp3.hpp
// /////////////////////////////////////////////////////////////////////////////
