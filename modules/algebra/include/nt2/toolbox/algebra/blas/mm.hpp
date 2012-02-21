/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_ALGEBRA_BLAS_MM_HPP_INCLUDED
#define NT2_TOOLBOX_ALGEBRA_BLAS_MM_B_HPP_INCLUDED
#include <nt2/toolbox/algebra/blas/blas3.hpp>

namespace nt2
{
  namespace details
  {
    //all mm call: gemm symm hemm trmm for data types float, double and related complex 
    // gemm C <- al*op(A)*op(B)+be*C (op t h n according ta and tb)
    // symm and hemm C <- al*A*B +be+C or   al*B*A +be+C (A sy or he and according side) 
    // trmm B <- al*op(A)*B or  B <- al*B*op(A) (op t h n) A triangular accessed according uplo
    
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

#define NT2_MM(T, PREFIX)                                               \
    inline void trmm(                                                   \
                     const char *side, const char *uplo,                \
                     const char *transa,                                \
                     const long int *m,                                 \
                     const long int *n,                                 \
                     const T *al, const T *a, const long int *lda,      \
                     T *b, const long int *ldb)                         \
    {
      BOOST_PP_CAT(PREFIX,BOOST_PP_CAT(trmm,_))(side,uplo,transa,m,n,al,a,lda,b,ldb); \
    }                                                                   \
    /**/

    NT2_MM(double, d) 
    NT2_MM(float,  s) 
    NT2_MM(std::complex<double>, z)
    NT2_MM(std::complex<float>, c)

#undef NT2_MM          
  }
}


#endif

