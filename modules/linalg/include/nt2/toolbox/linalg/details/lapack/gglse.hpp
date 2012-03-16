/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2009-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_DETAILS_LAPACK_GGLSE_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_DETAILS_LAPACK_GGLSE_HPP_INCLUDED
#include <nt2/toolbox/linalg/details/utility/f77_wrapper.hpp>
#include <nt2/toolbox/linalg/details/lapack/workspace.hpp>

extern "C"
{
  void NT2_F77NAME(cgglse)(const la_int* m, const la_int* n,
                           const la_int* p, la_complex* a, const la_int* lda, la_complex* b,
                           const la_int* ldb, la_complex* c, const la_complex* d, la_complex* x,
                           la_complex* work, const la_int* lwork, la_int* info);
  void NT2_F77NAME(zgglse)(const la_int* m, const la_int* n,
                           const la_int* p, la_complex* a, const la_int* lda, la_complex* b,
                           const la_int* ldb, la_complex* c, const la_complex* d, la_complex* x,
                           la_complex* work, const la_int* lwork, la_int* info);
  void NT2_F77NAME(dgglse)(const la_int* m, const la_int* n,
                           const la_int* p, double* a, const la_int* lda, double* b,
                           const la_int* ldb, double* c, const double* d, double* x,
                           double* work, const la_int* lwork, la_int* info);
  void NT2_F77NAME(sgglse)(const la_int* m, const la_int* n,
                           const la_int* p, float* a, const la_int* lda, float* b,
                           const la_int* ldb, float* c, const float* d, float* x,
                           float* work, const la_int* lwork, la_int* info);
}

namespace nt2
{
  namespace details
  {
#define NT2_GGLSE(NAME, T)                                          \
    inline void gglse( const la_int* m,                             \
                       const la_int* n,                             \
                       const la_int* p,                             \
                       T* a, const la_int* lda,                     \
                       T* b, const la_int* ldb,                     \
                       T* c,                                        \
                       const T* d,                                  \
                       T* x,                                        \
                       la_int* info,                                \
                       nt2::details::workspace<T> & w)              \
    {                                                               \
      NT2_F77NAME(NAME)(m, n, p, a, lda, b, ldb,                    \
                        c, d, x, w.getw(), w.query(), info);        \
      w.resizew(w.neededsize());                                    \
      NT2_F77NAME(NAME)(m, n, p, a, lda, b, ldb,                    \
                        c, d, x, w.getw(), &w.neededsize(), info);  \
    }                                                               \
    inline void gglse( const la_int* m,                             \
                       const la_int* n,                             \
                       const la_int* p,                             \
                       T* a, const la_int* lda,                     \
                       T* b, const la_int* ldb,                     \
                       T* c,                                        \
                       const T* d,                                  \
                       T* x,                                        \
                       la_int* info)                                \
    {                                                               \
      nt2::details::workspace<T> w;                                 \
      gglse(m, n, p, a, lda, b, ldb, c, d, x, info, w);             \
    }                                                               \
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
