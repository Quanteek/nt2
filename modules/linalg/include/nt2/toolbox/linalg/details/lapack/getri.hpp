/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_DETAILS_LAPACK_GETRI_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_DETAILS_LAPACK_GETRI_HPP_INCLUDED
#include <nt2/toolbox/linalg/details/utility/f77_wrapper.hpp>
#include <nt2/toolbox/linalg/details/lapack/workspace.hpp>
// trgetri

namespace nt2
{
  namespace details
  {
    extern "C"
    {
#define NT2_COMPLEX void
  void NT2_F77NAME(cgetri)(const long int* n, NT2_COMPLEX* a, const long int* lda,
                       const long int* ipiv, NT2_COMPLEX* work, const long int* lwork,
                       long int* info);
  void NT2_F77NAME(dgetri)(const long int* n, double* a, const long int* lda,
                       const long int* ipiv, double* work, const long int* lwork,
                       long int* info);
  void NT2_F77NAME(sgetri)(const long int* n, float* a, const long int* lda,
                       const long int* ipiv, float* work, const long int* lwork,
                       long int* info);
  void NT2_F77NAME(zgetri)(const long int* n, NT2_COMPLEX* a, const long int* lda,
                       const long int* ipiv, NT2_COMPLEX* work, const long int* lwork,
                       long int* info);
      
#undef NT2_COMPLEX
    }
    
#define NT2_GETRI(NAME, T)                      \
    inline void getri(const long int* n,        \
                      T* a,                     \
                      const long int* lda,      \
                      const long int* ipiv,     \
                      long int* info,           \
                      workspace<T> & w)         \
    {                                           \
      NT2_F77NAME(NAME)(n,a,lda,ipiv,           \
                        w.getw(),w.query(),     \
                        info);                  \
      w.resizew(w.neededsize());                \
      NT2_F77NAME(NAME)(n,a,lda,                \
                        ipiv,w.getw(),          \
                        &w.neededsize(),        \
                        info);                  \
    }                                           \
    inline void getri(const long int* n,        \
                      T* a,                     \
                      const long int* lda,      \
                      const long int* ipiv,     \
                      long int* info)           \
    {                                           \
      workspace<T> w;                           \
      getri(n, a, lda, ipiv, info, w);          \
    }                                           \


    NT2_GETRI(sgetri, float)
    NT2_GETRI(dgetri, double)
    NT2_GETRI(cgetri, std::complex<float>)
    NT2_GETRI(zgetri, std::complex<double>)
      
#undef NT2_GETRI
  }
}

#endif

// /////////////////////////////////////////////////////////////////////////////
// End of getri.hpp
// /////////////////////////////////////////////////////////////////////////////
