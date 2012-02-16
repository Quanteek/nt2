//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II       
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI       
//                                                                              
//          Distributed under the Boost Software License, Version 1.0.          
//                 See accompanying file LICENSE.txt or copy at                 
//                     http://www.boost.org/LICENSE_1_0.txt                     
//==============================================================================
#ifndef NT2_TOOLBOX_ALGEBRA_FUNCTIONS_BLAS_GENERAL_GEMV_HPP_INCLUDED
#define NT2_TOOLBOX_ALGEBRA_FUNCTIONS_BLAS_GENERAL_GEMV_HPP_INCLUDED

#include <nt2/table.hpp>
#include <nt2/core/container/category.hpp>
#include <nt2/include/functions/size.hpp>
#include <nt2/include/functions/extent.hpp>
#include <nt2/include/functions/numel.hpp>
#include <nt2/toolbox/algebra/blas/blas2.hpp>
#include <nt2/toolbox/algebra/details/padding.hpp>
#include <nt2/include/constants/zero.hpp>
#include <nt2/include/constants/one.hpp>
#include <nt2/include/functions/isvector.hpp>
#include <nt2/include/functions/iscolumn.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/assert.hpp>

namespace nt2
{

  template<char T0, char T1>
  struct gemv_status
  {
    static const char tA = T0;
    static const char tx = T1;
  };

  namespace ext
  {
    
#define NT2_GEMV(T, PREFIX)                                             \
    inline void gemv(const char *ta, const long int *m,                 \
                     const long int *n,                                 \
                     const T *al,                                       \
                     const T *a, const long int *lda,                   \
                     const T *x, const long int *incx,                  \
                     const T *be,                                       \
                     T *y, const long int *incy)                        \
    {                                                                   \
      BOOST_PP_CAT(PREFIX,BOOST_PP_CAT(gemv,_))(ta,m,n,al,a,lda,x,incx,be,y,incy);       \
    }                                                                   \
        
    NT2_GEMV(double, d)
    NT2_GEMV(float,  s)
    NT2_GEMV(std::complex<double>, z)
    NT2_GEMV(std::complex<float>, c)

#undef NT2_GEMV
      
    /////////////////////////////////////////////////////////////////////////////
    // Implementation when table_ are floating_
    /////////////////////////////////////////////////////////////////////////////
      NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::gemv_, tag::cpu_
                                , (A0)(S0)(A1)(S1)(A2)(S2)(A3)(A4)(A5)
                                , (unspecified_ < A5 > )
                                  ((table_< floating_<A0>, S0 > ))
                                  ((table_< floating_<A1>, S1 > ))
                                  ((table_< floating_<A2>, S2 > ))
                                  (scalar_ < arithmetic_<A3 > > )
                                  (scalar_ < arithmetic_<A4 > > )
                                  )
    {
      typedef void result_type;
      BOOST_FORCEINLINE result_type operator()( A5 const& 
                                              , A0 const& a, A1 const& x, A2& y
                                              , A3 const& al, A4 const& be
                                              )
      {
        typedef typename A0::value_type value_type; 
        const char transa = A5::tA; 
        const char transx = A5::tx; 

        BOOST_ASSERT_MSG( (is_vector(x) && is_column(y)), 
                          "Matrix-vector product y = al*A*x+be*y (gemv) must be called with vector x and column vector y");

        const long int incx = 1; //(transx == 'N') ? 1 : nt2::details::padding(x)/sizeof(value_type); 
        const long int incy = 1; 
        const long int m   = nt2::extent(a)[transa=='N'?0:1]; 
        const long int n   = nt2::extent(a)[transa=='N'?1:0];
        BOOST_ASSERT_MSG( (n == nt2::numel(x)), 
                          "In matrix-vector product y = al*A*x+be*y (gemv) the number of columns \
                          (lines if transpose) of A must be equal to the numel of x ");
        BOOST_ASSERT_MSG( (nt2::numel(y) == m),  
                           "In Matrix-vector product y = al*A*x+be*y (gemv) the number of lines \
                           (columns if transpose) of A must be equal to the number of lines of y");
        
        const value_type alpha = al; 
        const long int lda = nt2::details::padding(a);
        const value_type beta  = be;
        gemv(&transa,&m,&n,&alpha,a.begin(),&lda,x.begin(),&incx,&beta,y.begin(),&incy);
      }
    };
    
  }
  
  template < class A5, class A0, class A1, class A2, class A3>
  inline void gemv(A5 const& a5, A0 const& a0, A1 const& a1, A2& a2, A3 const& a3)
  {
    typedef typename A0::value_type value_type; 
    gemv(a5, a0, a1, a2, a3, Zero<value_type>());
  }
  
  template < class A5, class A0, class A1, class A2>
  inline void gemv(A5 const& a5, A0 const& a0, A1 const& a1, A2& a2)
  {
    typedef typename A0::value_type value_type; 
    gemv(a5, a0, a1, a2, One<value_type>(), Zero<value_type>());
  }
  template < char transa, char transx, class A0, class A1, class A2, class A3, class A4>
  inline void gemv(A0 const& a0, A1 const& a1, A2& a2, A3 const& a3, A4 const& a4)
  {
    gemv(gemv_status<transa, transx>(), a0, a1, a2, a3, a4);
  }
  
  template < char transa, char transx, class A0, class A1, class A2, class A3>
  inline void gemv(A0 const& a0, A1 const& a1, A2& a2, A3 const& a3)
  {
    typedef typename A0::value_type value_type; 
    gemv(gemv_status<transa, transx>(), a0, a1, a2, a3, Zero<value_type>());
  }
  
  template <char transa, char transx, class A0, class A1, class A2>
  inline void gemv(A0 const& a0, A1 const& a1, A2& a2)
  {
    typedef typename A0::value_type value_type; 
    gemv(gemv_status<transa, transx>(), a0, a1, a2, One<value_type>(), Zero<value_type>());
  }
  
}

#endif
