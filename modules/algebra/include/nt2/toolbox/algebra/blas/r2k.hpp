/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_ALGEBRA_BLAS_R2K_HPP_INCLUDED
#define NT2_TOOLBOX_ALGEBRA_BLAS_R2K_B_HPP_INCLUDED
#include <nt2/toolbox/algebra/blas/blas3.hpp>

namespace nt2
{
  namespace details
  {
    // tA is the transpose of A
    // hA is the transconjugate of A
    // all r2k call: ger2k syr2k her2k for data types float, double and related complex 
    // syr2k  C <- al*tA*B+al*B*tA+be*C or  C <- al*tA*B+A*tB+be*C (trans N, T)
    // her2k  C <- al*hA*B+al*B*hA+be*C or  C <- al*hA*B+A*hB+be*C (trans N, C)
    // The produced C is symetric or hermitian with uplo part filled

#define NT2_R2K(T, PREFIX)                                              \
    inline void syr2k(                                                  \
                     const char *uplo,                                  \
                     const long int *n,                                 \
                     const T *al,                                       \
                     const T *a, const long int *lda,                   \
                     const T *b, const long int *ldb,                   \
                     const T *be, T *c,                                 \
                     const long int *ldc)                               \
    {                                                                   \
      BOOST_PP_CAT(PREFIX,BOOST_PP_CAT(syr2k,_))(uplo,n,al,a,lda,be,c,ldc); \
    }                                                                   \
          /**/
      
    NT2_R2K(double, d) 
    NT2_R2K(float,  s) 
    NT2_R2K(std::complex<double>, z) 
    NT2_R2K(std::complex<float>, c) 

#undef NT2_R2K

#define NT2_R2K(T, PREFIX)                                              \
      inline void her2k(                                                \
                     const char *uplo,                                  \
                     const long int *n,                                 \
                     const T *al,                                       \
                     const T *a, const long int *lda,                   \
                     const T *b, const long int *ldb,                   \
                     const T *be, T *c,                                 \
                     const long int *ldc)                               \
    {                                                                   \
      BOOST_PP_CAT(PREFIX,BOOST_PP_CAT(her2k,_))(uplo,n,al,a,lda,be,c,ldc); \
    }                                                                   \
    /**/

    NT2_R2K(std::complex<double>, z)
    NT2_R2K(std::complex<float>, c)

#undef NT2_R2K
  }
}


#endif

