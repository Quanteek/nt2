/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_DETAILS_BLAS_SV_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_DETAILS_BLAS_SV_HPP_INCLUDED
#include <nt2/toolbox/linalg/details/blas/blas2.hpp>

// all sv call:  trsv  tbsv etc.
// for data types float, double and related complex 
// _TRSV ( UPLO, TRANS, DIAG,    N,                A, LDA, X, INCX )                 S, D, C, Z
// _TBSV ( UPLO, TRANS, DIAG,    N, K,             A, LDA, X, INCX )                 S, D, C, Z
// _TPSV ( UPLO, TRANS, DIAG,    N,                AP,     X, INCX )                 S, D, C, Z



namespace nt2
{
  namespace details
  {

#define NT2_SV(T, PREFIX)                                               \
    inline void trsv(const char *uplo, const char *trans,               \
                     const char *diag,                                  \
                     const long int *n,                                 \
                     const T *a, const long int *lda                    \
                     T *x, const long int *incx)                        \
    {                                                                   \
      BOOST_PP_CAT(PREFIX,BOOST_PP_CAT(tpsv,_))                         \
        (uplo,trans,diag,n,a,lda,x,incx);                               \
    }                                                                   \
        
      NT2_SV(double, d)
      NT2_SV(float,  s)
      NT2_SV(std::complex<double>, z)
      NT2_SV(std::complex<float>, c)
      
#undef NT2_SV
      
#define NT2_SV(T, PREFIX)                                               \
    inline void tbsv(const char *uplo, const char *trans,               \
                     const char *diag,                                  \
                     const long int *n,  const long int *k,             \
                     const T *a,                                        \
                     const long int *lda,                               \
                     T *x, const long int *incx)                        \
    {                                                                   \
      BOOST_PP_CAT(PREFIX,BOOST_PP_CAT(trsv,_))                         \
        (uplo,trans,diag,n,k, a,lda,x,incx);                            \
    }                                                                   \
        
      NT2_SV(double, d)
      NT2_SV(float,  s)
      NT2_SV(std::complex<double>, z)
      NT2_SV(std::complex<float>, c)
      
#undef NT2_SV
        
#define NT2_SV(T, PREFIX)                                               \
    inline void tpsv(const char *uplo, const char *trans,               \
                     const char *diag,                                  \
                     const long int *n,                                 \
                     const T *ap,                                       \
                     T *x, const long int *incx)                        \
    {                                                                   \
      BOOST_PP_CAT(PREFIX,BOOST_PP_CAT(tpsv,_))                         \
        (uplo,trans,diag,n,ap,x,incx);                               \
    }                                                                   \
        
      NT2_SV(double, d)
      NT2_SV(float,  s)
      NT2_SV(std::complex<double>, z)
      NT2_SV(std::complex<float>, c)
      
#undef NT2_SV


        
  }  
}


#endif

// /////////////////////////////////////////////////////////////////////////////
// End of sv.hpp
// /////////////////////////////////////////////////////////////////////////////
