/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_DETAILS_BLAS_BLAS3_B_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_DETAILS_BLAS_BLAS3_B_HPP_INCLUDED
#include <nt2/toolbox/linalg/details/blas/blas3.hpp>

namespace nt2
{
  namespace details
  {
    //all mm call: gemm symm hemm for data types float, double and related complex 
    
#define NT2_MM(T, PREFIX)                                               \
    inline void gemm(const char *ta, const char *tb, const long int *m, \
                     const long int *n, const long int *k,              \
                     const T *al, const T *a,                           \
                     const long int *lda, const T *b,                   \
                     const long int *ldb, const T *be, T *c,            \
                     const long int *ldc)                               \
    {                                                                   \
      BOOST_PP_CAT(PREFIX,BOOST_PP_CAT(gemm,_))(ta,tb,m,n,k,al,a,lda,b,ldb,be,c,ldc); \
    }                                                                   \
    /**/

    NT2_MM(double, d) 
    NT2_MM(float,  s) 
    NT2_MM(std::complex<double>, z) 
    NT2_MM(std::complex<float>, c) 

#undef NT2_MM
      
#define NT2_MM(T, PREFIX)                                               \
    inline void symm(                                                   \
                   const char *side, const char *uplo,                  \
                   const long int *m,                                   \
                   const long int *n,                                   \
                   const T *al, const T *a,                             \
                   const long int *lda, const T *b,                     \
                   const long int *ldb, const T *be, T *c,              \
                   const long int *ldc)                                 \
    {                                                                   \
      BOOST_PP_CAT(PREFIX,BOOST_PP_CAT(symm,_))(side,uplo,m,n,al,a,lda,b,ldb,be,c,ldc); \
    }                                                                   \
          /**/
      
    NT2_MM(double, d) 
    NT2_MM(float,  s) 
    NT2_MM(std::complex<double>, z) 
    NT2_MM(std::complex<float>, c) 

#undef NT2_MM

#define NT2_MM(T, PREFIX)                                               \
    inline void hemm(                                                   \
                     const char *side, const char *uplo,                \
                     const long int *m,                                 \
                     const long int *n,                                 \
                     const T *al, const T *a,                           \
                     const long int *lda, const T *b,                   \
                     const long int *ldb, const T *be, T *c,            \
                     const long int *ldc)                               \
    {                                                                   \
      BOOST_PP_CAT(PREFIX,BOOST_PP_CAT(hemm,_))(side,uplo,m,n,al,a,lda,b,ldb,be,c,ldc); \
    }                                                                   \
    /**/

    NT2_MM(std::complex<double>, z)
    NT2_MM(std::complex<float>, c)

#undef NT2_MM
  }
}


#endif

