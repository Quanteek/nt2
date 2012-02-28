/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2009-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_DETAILS_LAPACK_LSE_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_DETAILS_LAPACK_LSE_HPP_INCLUDED
#include <nt2/toolbox/linalg/details/utility/f77_wrapper.hpp>
#include <nt2/toolbox/linalg/details/lapack/workspace.hpp>

namespace nt2
{
  namespace details
  {
    extern "C"
    {
      #define NT2_COMPLEX void
      void NT2_F77NAME(cgglse)(const long int* m, const long int* n,
                           const long int* p, NT2_COMPLEX* a, const long int* lda, NT2_COMPLEX* b,
                           const long int* ldb, NT2_COMPLEX* c, const NT2_COMPLEX* d, NT2_COMPLEX* x,
                           NT2_COMPLEX* work, const long int* lwork, long int* info);
      void NT2_F77NAME(zgglse)(const long int* m, const long int* n,
                           const long int* p, NT2_COMPLEX* a, const long int* lda, NT2_COMPLEX* b,
                           const long int* ldb, NT2_COMPLEX* c, const NT2_COMPLEX* d, NT2_COMPLEX* x,
                           NT2_COMPLEX* work, const long int* lwork, long int* info);
      void NT2_F77NAME(dgglse)(const long int* m, const long int* n,
                           const long int* p, double* a, const long int* lda, double* b,
                           const long int* ldb, double* c, const double* d, double* x,
                           double* work, const long int* lwork, long int* info);
      void NT2_F77NAME(sgglse)(const long int* m, const long int* n,
                           const long int* p, float* a, const long int* lda, float* b,
                           const long int* ldb, float* c, const float* d, float* x,
                           float* work, const long int* lwork, long int* info);
      #undef NT2_COMPLEX
    }
    
#define NT2_GGLSE(NAME, T)                                              \
    inline void gglse( const long int* m,                               \
                       const long int* n,                               \
                       const long int* p,                               \
                       T* a, const long int* lda,                       \
                       T* b, const long int* ldb,                       \
                       T* c,                                            \
                       const T* d,                                      \
                       T* x,                                            \
                       long int* info,                                  \
                       ::nt2::details::workspace<T> & w)                \
    {                                                                   \
      NT2_F77NAME(NAME)(m, n, p, a, lda, b, ldb,                        \
                        c, d, x, w.getw(), w.query(), info);            \
      w.resizew(w.neededsize());                                        \
      NT2_F77NAME(NAME)(m, n, p, a, lda, b, ldb,                        \
                        c, d, x, w.getw(), &w.neededsize(), info);      \
    }                                                                   \
    inline void gglse( const long int* m,                               \
                       const long int* n,                               \
                       const long int* p,                               \
                       T* a, const long int* lda,                       \
                       T* b, const long int* ldb,                       \
                       T* c,                                            \
                       const T* d,                                      \
                       T* x,                                            \
                       long int* info)                                  \
    {                                                                   \
      ::nt2::details::workspace<T> w;                                   \
      gglse(m, n, p, a, lda, b, ldb, c, d, x, info, w);                 \
    }                                                                   \
        /**/
        
    NT2_GGLSE(sgglse, float)
    NT2_GGLSE(dgglse, double)
    NT2_GGLSE(cgglse, std::complex<float>)
    NT2_GGLSE(zgglse, std::complex<double>)

#undef NT2_GGLSE
    
  }
}


#endif

// /////////////////////////////////////////////////////////////////////////////
// End of l_lse.hpp<3>
// /////////////////////////////////////////////////////////////////////////////
