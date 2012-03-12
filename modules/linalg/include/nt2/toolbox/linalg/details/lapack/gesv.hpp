/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_DETAILS_LAPACK_GESV_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_DETAILS_LAPACK_GESV_HPP_INCLUDED
#include <nt2/toolbox/linalg/details/utility/f77_wrapper.hpp>
namespace nt2
{
  namespace details
  {
    extern "C"
    {
#define NT2_COMPLEX void
  void NT2_F77NAME(cgesv)(const long int* n, const long int* nrhs,
                          NT2_COMPLEX* a, const long int* lda, long int* ipiv,
                          NT2_COMPLEX* b, const long int* ldb, long int* info);
  void NT2_F77NAME(dgesv)(const long int* n, const long int* nrhs,
                          double* a, const long int* lda, long int* ipiv,
                          double* b, const long int* ldb, long int* info);
  void NT2_F77NAME(sgesv)(const long int* n, const long int* nrhs,
                          float* a, const long int* lda, long int* ipiv,
                          float* b, const long int* ldb, long int* info);
  void NT2_F77NAME(zgesv)(const long int* n, const long int* nrhs,
                          NT2_COMPLEX* a, const long int* lda, long int* ipiv,
                          NT2_COMPLEX* b, const long int* ldb, long int* info);
#undef NT2_COMPLEX
    }
    
#define NT2_GESV(NAME, T)                       \
    inline void gesv(const long int* n,         \
                     const long int* nrhs,      \
                     T* a,                      \
                     const long int* lda,       \
                     long int* ipiv,            \
                     T* b,                      \
                     const long int* ldb,       \
                     long int* info,            \
                     workspace<T> & w)          \
    {                                           \
      NT2_F77NAME( NAME )(n, nrhs, a, lda,      \
                          ipiv, b, ldb, info);  \
    }                                           \
    inline void gesv(const long int* n,         \
                     const long int* nrhs,      \
                     T* a,                      \
                     const long int* lda,       \
                     long int* ipiv,            \
                     T* b,                      \
                     const long int* ldb,       \
                     long int* info)            \
    {                                           \
      workspace<T> w;                           \
      gesv(n, nrhs, a, lda,                     \
           ipiv, b, ldb, info, w);              \
    }                                           \

        
    NT2_GESV(cgesv, float)
    NT2_GESV(zgesv, double)
    NT2_GESV(cgesv, std::complex<float>)
    NT2_GESV(zgesv, std::complex<double>)
      
#undef NT2_GESV
  }
}

#endif

// /////////////////////////////////////////////////////////////////////////////
// End of svx.hpp
// /////////////////////////////////////////////////////////////////////////////
