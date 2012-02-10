//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II       
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI       
//                                                                              
//          Distributed under the Boost Software License, Version 1.0.          
//                 See accompanying file LICENSE.txt or copy at                 
//                     http://www.boost.org/LICENSE_1_0.txt                     
//==============================================================================
#ifndef NT2_TOOLBOX_ALGEBRA_FUNCTIONS_BLAS_GEMV_HPP_INCLUDED
#define NT2_TOOLBOX_ALGEBRA_FUNCTIONS_BLAS_GEMV_HPP_INCLUDED

#include <nt2/table.hpp>
#include <nt2/core/container/category.hpp>
#include <nt2/include/functions/size.hpp>
#include <nt2/toolbox/algebra/blas/blas2.hpp>
#include <boost/simd/sdk/memory/align_on.hpp>

#define F77NAME(x) x##_

namespace nt2 { namespace ext
{
/////////////////////////////////////////////////////////////////////////////
// Implementation when type A0 is table_<double_>
/////////////////////////////////////////////////////////////////////////////
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::gemv_, tag::cpu_
                            , (A0)(S0)(A1)(S1)(A2)(S2)
                            , ((table_< double_<A0>, S0 > ))
                              ((table_< double_<A1>, S1 > ))
                              ((table_< double_<A2>, S2 > ))
                            )
  {
    typedef void result_type;
    typedef typename A1::parent::lead_t lead_t_a1;

    BOOST_FORCEINLINE result_type operator()(A0& a0, A1 const& a1, A2 const& a2)
    {
      const char transa  = 'N';
      const long int m   = nt2::size(a1)(1);
      const long int n   = nt2::size(a1)(2);
      const double alpha = 1.0; 
      const long int lda = boost::simd::memory::align_on(m, lead_t_a1::value);
      const long int x   = 1;
      const double beta  = 0.0; 
      const long int y   = 1;
      F77NAME(dgemv)(&transa,&m,&n,&alpha,a1.begin(),&lda,a2.begin(),&x,&beta,a0.begin(),&y);
    }
  };

/////////////////////////////////////////////////////////////////////////////
// Implementation when type A0 is table_<float_>
/////////////////////////////////////////////////////////////////////////////
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::gemv_, tag::cpu_
                            , (A0)(S0)(A1)(S1)(A2)(S2)
                            , ((table_< single_<A0>, S0 > ))
                              ((table_< single_<A1>, S1 > ))
                              ((table_< single_<A2>, S2 > ))
                            )
  {
    typedef void result_type;
    typedef typename A1::parent::lead_t lead_t_a1;
  
    BOOST_FORCEINLINE result_type operator()(A0& a0, A1 const& a1, A2 const& a2)
    {
      const char transa  = 'N';
      const long int m   = nt2::size(a1)(1);
      const long int n   = nt2::size(a1)(2);
      const float alpha  = 1.0; 
      const long int lda = boost::simd::memory::align_on(m, lead_t_a1::value);
      const long int x   = 1;
      const float beta   = 0.0; 
      const long int y   = 1;
      F77NAME(sgemv)(&transa,&m,&n,&alpha,a1.begin(),&lda,a2.begin(),&x,&beta,a0.begin(),&y);
    }
  };

} }

#undef F77NAME

#endif
