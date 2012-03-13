/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_DETAILS_LAPACK_GELS_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_DETAILS_LAPACK_GELS_HPP_INCLUDED
#include <nt2/toolbox/linalg/details/utility/f77_wrapper.hpp>
#include <nt2/toolbox/linalg/details/lapack/workspace.hpp>
namespace nt2
{
  namespace details
  {
    extern "C"
    {
#define NT2_COMPLEX void
      void  NT2_F77NAME(cgels)(const char* trans, const long int* m, const long int* n, const long int* nrhs,
                               const NT2_COMPLEX* a, const long int* lda, NT2_COMPLEX* b, const long int* ldb,
                               NT2_COMPLEX* work, const long int* lwork, long int* info);
      void  NT2_F77NAME(dgels)(const char* trans, const long int* m, const long int* n, const long int* nrhs,
                               const double* a, const long int* lda, double* b, const long int* ldb,
                               double* work, const long int* lwork, long int* info);
      void  NT2_F77NAME(sgels)(const char* trans, const long int* m, const long int* n, const long int* nrhs,
                               const float* a, const long int* lda, float* b, const long int* ldb,
                               float* work, const long int* lwork, long int* info);
      void  NT2_F77NAME(zgels)(const char* trans, const long int* m, const long int* n, const long int* nrhs,
                               const NT2_COMPLEX* a, const long int* lda, NT2_COMPLEX* b, const long int* ldb,
                               NT2_COMPLEX* work, const long int* lwork, long int* info);
#undef NT2_COMPLEX
    }
    
#define NT2_GELS(NAME, T)                       \
    inline void gels(const char* trans,         \
                     const long int* m,         \
                     const long int* n,         \
                     const long int* nrhs,      \
                     const T* a,                \
                     const long int* lda,       \
                     T* b,                      \
                     const long int* ldb,       \
                     long int* info,            \
                     workspace<T> & w)          \
    {                                           \
      NT2_F77NAME( NAME )(trans, m, n, nrhs,    \
                      a, lda, b, ldb,           \
                      w.getw(), w.query(),      \
                      info);                    \
      w.resizew(w.neededsize());                \
      NT2_F77NAME( NAME )(trans, m, n, nrhs,    \
                      a, lda, b, ldb,           \
                      w.getw(),&w.neededsize(), \
                      info);                    \
    }                                           \
    inline void gels(const char* trans,         \
                     const long int* m,         \
                     const long int* n,         \
                     const long int* nrhs,      \
                     const T* a,                \
                     const long int* lda,       \
                     T* b,                      \
                     const long int* ldb,       \
                     long int* info)            \
    {                                           \
      workspace<T> w;                           \
      gels(trans, m, n, nrhs,                   \
           a, lda, b, ldb,                      \
           info, w);                            \
    }                                           \
        
    NT2_GELS(cgels, float)
    NT2_GELS(zgels, double)
    NT2_GELS(cgels, std::complex<float>)
    NT2_GELS(zgels, std::complex<double>)
      
#undef NT2_GELS
  }
}

#endif

