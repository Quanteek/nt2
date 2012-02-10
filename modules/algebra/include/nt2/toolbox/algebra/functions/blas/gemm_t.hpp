//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II       
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI       
//                                                                              
//          Distributed under the Boost Software License, Version 1.0.          
//                 See accompanying file LICENSE.txt or copy at                 
//                     http://www.boost.org/LICENSE_1_0.txt                     
//==============================================================================
#ifndef NT2_TOOLBOX_ALGEBRA_FUNCTIONS_BLAS_GEMM_T_HPP_INCLUDED
#define NT2_TOOLBOX_ALGEBRA_FUNCTIONS_BLAS_GEMM_T_HPP_INCLUDED

#include <nt2/table.hpp>
#include <nt2/core/container/category.hpp>
#include <nt2/include/functions/size.hpp>
#include <nt2/toolbox/algebra/blas/blas3.hpp>
#include <boost/simd/sdk/memory/align_on.hpp>

#define F77NAME(x) x##_

namespace nt2 { namespace ext
{
/////////////////////////////////////////////////////////////////////////////
// Implementation when type A0 is table_<double_>
/////////////////////////////////////////////////////////////////////////////
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::gemm_t_, tag::cpu_
                            , (A0)(S0)(A1)(S1)(A2)(S2)
                            , ((table_< double_<A0>, S0 > ))
                              ((table_< double_<A1>, S1 > ))
                              ((table_< double_<A2>, S2 > ))
                            )
  {
    typedef void result_type;
    typedef typename A1::parent::lead_t lead_t_a1;
    typedef typename A2::parent::lead_t lead_t_a2;

    BOOST_FORCEINLINE result_type operator()(A0& a0, A1 const& a1, A2 const& a2)
    {
      const char transa = 'N';
      const char transb = 'T';
      const long int m = nt2::size(a1)(1);
      const long int n = nt2::size(a2)(1);
      const long int k = nt2::size(a1)(2);
      const double alpha = 1.0; 
      const long int lda = boost::simd::memory::align_on(m, lead_t_a1::value);
      const long int ldb = boost::simd::memory::align_on(n, lead_t_a2::value);
      const double beta = 0.0; 
      const long int ldc = boost::simd::memory::align_on(m, lead_t_a1::value);
      F77NAME(dgemm)(&transa,&transb,&m,&n,&k,&alpha,a1.begin(),&lda,a2.begin(),&ldb,&beta,a0.begin(),&ldc);
    }
  };

/////////////////////////////////////////////////////////////////////////////
// Implementation when type A0 is table_<float_>
/////////////////////////////////////////////////////////////////////////////
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::gemm_t_, tag::cpu_
                            , (A0)(S0)(A1)(S1)(A2)(S2)
                            , ((table_< single_<A0>, S0 > ))
                              ((table_< single_<A1>, S1 > ))
                              ((table_< single_<A2>, S2 > ))
                            )
  {
    typedef void result_type;
    typedef typename A1::parent::lead_t lead_t_a1;
    typedef typename A2::parent::lead_t lead_t_a2;
    
    BOOST_FORCEINLINE result_type operator()(A0& a0, A1 const& a1, A2 const& a2)
    {
      const char transa = 'N';
      const char transb = 'T';
      const long int m = nt2::size(a1)(1);
      const long int n = nt2::size(a2)(1);
      const long int k = nt2::size(a1)(2);
      const float alpha = 1.0; 
      const long int lda = boost::simd::memory::align_on(m, lead_t_a1::value);
      const long int ldb = boost::simd::memory::align_on(n, lead_t_a2::value);
      const float beta = 0.0; 
      const long int ldc = boost::simd::memory::align_on(m, lead_t_a1::value);
      F77NAME(sgemm)(&transa,&transb,&m,&n,&k,&alpha,a1.begin(),&lda,a2.begin(),&ldb,&beta,a0.begin(),&ldc);
    }
  };

} }

#undef F77NAME

#endif
