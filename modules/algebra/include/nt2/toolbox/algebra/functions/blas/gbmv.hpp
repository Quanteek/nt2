//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II       
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI       
//                                                                              
//          Distributed under the Boost Software License, Version 1.0.          
//                 See accompanying file LICENSE.txt or copy at                 
//                     http://www.boost.org/LICENSE_1_0.txt                     
//==============================================================================
#ifndef NT2_TOOLBOX_ALGEBRA_FUNCTIONS_BLAS_GENERAL_GBMV_HPP_INCLUDED
#define NT2_TOOLBOX_ALGEBRA_FUNCTIONS_BLAS_GENERAL_GBMV_HPP_INCLUDED

#include <nt2/table.hpp>
#include <nt2/core/container/category.hpp>
#include <nt2/include/functions/size.hpp>
#include <nt2/toolbox/algebra/blas/blas2.hpp>
#include <nt2/toolbox/algebra/details/padding.hpp>
#include <nt2/include/constants/zero.hpp>
#include <nt2/include/constants/one.hpp>
#include <nt2/include/functions/isvector.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/assert.hpp>
namespace nt2
{
  namespace ext
  {
    
#define NT2_GBMV(T, PREFIX)                                             \
    inline void gbmv(const char *ta, const long int *m,                 \
                     const long int *n,                                 \
                     const long int *kl,                                \
                     const long int *ku,                                \
                     const T *al,                                       \
                     const T *a, const long int *lda,                   \
                     const T *x, const long int *incx,                  \
                     const T *be,                                       \
                     T *y, const long int *incy)                        \
    {                                                                   \
      BOOST_PP_CAT(PREFIX,BOOST_PP_CAT(gbmv,_))(ta,m,n,kl,ku,al,a,lda,x,incx,be,y,incy);       \
    }                                                                   \
        
    NT2_GBMV(double, d)
    NT2_GBMV(float,  s)
    NT2_GBMV(std::complex<double>, z)
    NT2_GBMV(std::complex<float>, c)

#undef NT2_GBMV
      
    /////////////////////////////////////////////////////////////////////////////
    // Implementation when table_ are floating_
    /////////////////////////////////////////////////////////////////////////////
      NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::gbmv_, tag::cpu_
                                , (A0)(S0)(A1)(S1)(A2)(S2)(A3)(A4)(A5)
                                , (unspecified_ < A5 > )
                                  ((band_< floating_<A0>, S0 > ))
                                  ((band_< floating_<A1>, S1 > ))
                                  ((band_< floating_<A2>, S2 > ))
                                  (scalar_ < arithmetic_<A3 > > )
                                  (scalar_ < arithmetic_<A4 > > )
                                  )
    {
      typedef void result_type;
      BOOST_FORCEINLINE result_type operator()( A5 const& transa
                                              , A0 const& a, A1 const& x, A2& y
                                              , A3 const& al, A4 const& be
                                              )
      {
        typedef typename A0::value_type value_type; 
        BOOST_ASSERT_MSG( (is_vector(x) && is_vector(y)), 
                          "Matrix-vector product y = al*A*x+be*y (gbmv) must be called with 1D table for x and y");
        BOOST_ASSERT_MSG( (nt2::numel(x) == nt2::numel(y)),  
                           "Matrix-vector product y = al*A*x+be*y (gbmv) must be called with x and y of same size");

        const long int m   = nt2::size(a)(transa=='T'?2:1);
        const long int n   = nt2::size(a)(transa=='T'?1:2);

        BOOST_ASSERT_MSG( (n == nt2::numel(x)), 
                          "In matrix-vector product y = al*A*x+be*y (gemv) the inner number of columns \
                          (lines if transpose) of A must be equal to the numel of x ");

        const long int m   = nt2::size(a)(transa=='T'?2:1);
        const long int n   = nt2::size(a)(transa=='T'?1:2);
        const value_type alpha = al; 
        const long int lda = nt2::details::padding(a);
        const value_type beta  = be;
        const long int incx = 1; 
        const long int incy = 1;
        const long int kl = a.kl();
        const long int ku = a.ku(); 
        gbmv(&transa,&m,&n,&kl,&ku,&alpha,a.begin(),&lda,x.begin(),&incx,&beta,y.begin(),&incy);
      }
    };
    
  }

  template < class A5,  class A0,  class A1,  class A2,  class A3>
  inline void gbmv(A5 const& a5, A0 const& a0, A1 const& a1, A2& a2,A3 const& a3)
  {
    typedef typename A0::value_type value_type; 
    gbmv(a5, a0, a1, a2, a3, Zero<value_type>());
  }
  
  template < class A5,  class A0,  class A1,  class A2>
  inline void gbmv(A5 const& a5, A0 const& a0, A1 const& a1, A2& a2)
  {
    typedef typename A0::value_type value_type; 
    gbmv(a5, a0, a1, a2, One<value_type>(), Zero<value_type>());
  }
}

#endif
