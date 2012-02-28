/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_DETAILS_BLAS_RK_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_DETAILS_BLAS_RK_B_HPP_INCLUDED
#include <nt2/toolbox/linalg/details/blas/blas3.hpp>

namespace nt2
{
  namespace details
  {
    // tA is the transpose of A
    // hA is the transconjugate of A
    // all rk call: gerk syrk herk for data types float, double and related complex 
    // syrk  C <- al*A*tA+be*C or  C <- al*tA*A+be*C (trans N, T)
    // herk  C <- al*A*hA+be*C or  C <- al*A*hA+be*C (trans N, C)
    

#define NT2_RK(T, PREFIX)                                               \
    inline void syrk(                                                   \
                     const char *uplo,                                  \
                     const long int *n,                                 \
                     const T *al, const T *a,                           \
                     const long int *lda,                               \
                     const T *be, T *c,                                 \
                     const long int *ldc)                               \
    {                                                                   \
      BOOST_PP_CAT(PREFIX,BOOST_PP_CAT(syrk,_))(uplo,n,al,a,lda,be,c,ldc); \
    }                                                                   \
          /**/
      
    NT2_RK(double, d) 
    NT2_RK(float,  s) 
    NT2_RK(std::complex<double>, z) 
    NT2_RK(std::complex<float>, c) 

#undef NT2_RK

#define NT2_RK(T, PREFIX)                                               \
    inline void herk(                                                   \
                     const char *uplo,                                  \
                     const long int *n,                                 \
                     const T *al, const T *a,                           \
                     const long int *lda,                               \
                     const T *be, T *c,                                 \
                     const long int *ldc)                               \
    {                                                                   \
      BOOST_PP_CAT(PREFIX,BOOST_PP_CAT(herk,_))(uplo,n,al,a,lda,be,c,ldc); \
    }                                                                   \
    /**/

    NT2_RK(std::complex<double>, z)
    NT2_RK(std::complex<float>, c)

#undef NT2_RK
  }
}


#endif

