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
#include <nt2/toolbox/linalg/details/lapack/workspace.hpp>

extern "C"
{
  void NT2_F77NAME(cgesv)(const long int* n, const long int* nrhs,
                          la_complex* a, const long int* lda, long int* ipiv,
                          la_complex* b, const long int* ldb, long int* info);
  void NT2_F77NAME(dgesv)(const long int* n, const long int* nrhs,
                          double* a, const long int* lda, long int* ipiv,
                          double* b, const long int* ldb, long int* info);
  void NT2_F77NAME(sgesv)(const long int* n, const long int* nrhs,
                          float* a, const long int* lda, long int* ipiv,
                          float* b, const long int* ldb, long int* info);
  void NT2_F77NAME(zgesv)(const long int* n, const long int* nrhs,
                          la_complex* a, const long int* lda, long int* ipiv,
                          la_complex* b, const long int* ldb, long int* info);
}

namespace nt2
{
  namespace details
  {
    
#define NT2_GESV(NAME, T)                       \
    inline void gesv(const long int* n,         \
                     const long int* nrhs,      \
                     T* a,                      \
                     const long int* lda,       \
                     long int* ipiv,            \
                     T* b,                      \
                     const long int* ldb,       \
                     long int* info,            \
                     nt2::details::workspace<T> & )    \
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
      NT2_F77NAME( NAME )(n, nrhs, a, lda,      \
                          ipiv, b, ldb, info);  \
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
