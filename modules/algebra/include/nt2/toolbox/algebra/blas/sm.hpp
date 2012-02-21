/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_ALGEBRA_BLAS_SM_HPP_INCLUDED
#define NT2_TOOLBOX_ALGEBRA_BLAS_SM_B_HPP_INCLUDED
#include <nt2/toolbox/algebra/blas/blas3.hpp>

namespace nt2
{
  namespace details
  {
    //all sm call: trsm for data types float, double and related complex 

#define NT2_SM(T, PREFIX)                                               \
    inline void trsm(                                                   \
                     const char *side, const char *uplo,                \
                     const char *transa,                                \
                     const long int *m,                                 \
                     const long int *n,                                 \
                     const T *al, const T *a, const long int *lda,      \
                     T *b, const long int *ldb)                         \
    {
      BOOST_PP_CAT(PREFIX,BOOST_PP_CAT(trsm,_))(side,uplo,transa,m,n,al,a,lda,b,ldb); \
    }                                                                   \
    /**/

    NT2_SM(double, d) 
    NT2_SM(float,  s) 
    NT2_SM(std::complex<double>, z)
    NT2_SM(std::complex<float>, c)

#undef NT2_SM      
  }
}


#endif

