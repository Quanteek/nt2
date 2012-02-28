/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2009-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_DETAILS_LAPACK_LS_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_DETAILS_LAPACK_LS_HPP_INCLUDED
#include <nt2/toolbox/linalg/details/utility/f77_wrapper.hpp>
#include <nt2/toolbox/linalg/details/lapack/workspace.hpp>

namespace nt2
{
  namespace details
  {
    extern "C"
    {
      #define NT2_COMPLEX void
      void NT2_F77NAME(cgelsy)(const long int* m, const long int* n, const long int* nrhs,
                           NT2_COMPLEX* a, const long int* lda, NT2_COMPLEX* b,
                           const long int* ldb, long int* jpvt, const float* rcond,
                           long int* rank, NT2_COMPLEX* work, const long int* lwork, float* rwork, long int* info);
      void NT2_F77NAME(zgelsy)(const long int* m, const long int* n, const long int* nrhs,
                           NT2_COMPLEX* a, const long int* lda, NT2_COMPLEX* b,
                           const long int* ldb, long int* jpvt, const double* rcond,
                           long int* rank, NT2_COMPLEX* work, const long int* lwork, double* rwork, long int* info);
      void NT2_F77NAME(dgelsy)(const long int* m, const long int* n, const long int* nrhs,
                           double* a, const long int* lda, double* b,
                           const long int* ldb, long int* jpvt, const double* rcond,
                           long int* rank, double* work, const long int* lwork, long int* info);
      void NT2_F77NAME(sgelsy)(const long int* m, const long int* n, const long int* nrhs,
                           float* a, const long int* lda, float* b,
                           const long int* ldb, long int* jpvt, const float* rcond,
                           long int* rank, float* work, const long int* lwork, long int* info);
      #undef NT2_COMPLEX
    }
    
#define NT2_GELSY(NAME, T)                      \
    inline void gelsy(const long int* m,        \
                      const long int* n,        \
                      const long int* nrhs,     \
                      T* a,                     \
                      const long int* lda,      \
                      T* b,                     \
                      const long int* ldb,      \
                      long int* jpvt,           \
                      const T* rcond,           \
                      long int* rank,           \
                      long int* info,           \
                      nt2::details::workspace<T> & w)                   \
    {                                                                   \
      NT2_F77NAME( NAME )(m, n, nrhs, a, lda, b, ldb, jpvt,             \
                          rcond, rank, w.getw(), w.query(), info);      \
      w.resizew(w.neededsize());                                        \
      NT2_F77NAME( NAME )(m, n, nrhs, a, lda, b, ldb, jpvt,             \
                          rcond, rank, w.getw(), &w.neededsize(), info);\
    }                                                                   \
        inline void gelsy(const long int* m,                            \
                          const long int* n,                            \
                          const long int* nrhs,                         \
                          T* a,                                         \
                          const long int* lda,                          \
                          T* b,                                         \
                          const long int* ldb,                          \
                          long int* jpvt,                               \
                          const T* rcond,                               \
                          long int* rank,                               \
                          long int* info)                               \
        {                                                               \
          nt2::details::workspace<T> w;                                 \
          gelsy(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, info, w);\
        }                                                               \
        
    NT2_GELSY(sgelsy, float)
    NT2_GELSY(dgelsy, double)
    NT2_GELSY(cgelsy, std::complex<float>,  float)
    NT2_GELSY(zgelsy, std::complex<double>, double)

#undef NT2_GELSY
  }
}


#endif

// /////////////////////////////////////////////////////////////////////////////
// End of l_ls.hpp
// /////////////////////////////////////////////////////////////////////////////
