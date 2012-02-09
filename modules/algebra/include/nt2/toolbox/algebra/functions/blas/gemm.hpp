//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II       
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI       
//                                                                              
//          Distributed under the Boost Software License, Version 1.0.          
//                 See accompanying file LICENSE.txt or copy at                 
//                     http://www.boost.org/LICENSE_1_0.txt                     
//==============================================================================
#ifndef NT2_TOOLBOX_ALGEBRA_FUNCTIONS_BLAS_GEMM_HPP_INCLUDED
#define NT2_TOOLBOX_ALGEBRA_FUNCTIONS_BLAS_GEMM_HPP_INCLUDED

#include <nt2/table.hpp>
#include <nt2/core/functions/extent.hpp>
#include <nt2/core/container/category.hpp>
#include <nt2/include/functions/size.hpp>
#include <nt2/toolbox/algebra/blas/blas3.hpp>

namespace nt2 { namespace ext
{
/////////////////////////////////////////////////////////////////////////////
// Implementation when type A0 is table_<double_>
/////////////////////////////////////////////////////////////////////////////
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::gemm_, tag::cpu_
                            , (A0)(S0)(A1)(S1)(A2)(S2)
                            , ((table_< double_<A0>, S0 > ))
                              ((table_< double_<A1>, S1 > ))
                              ((table_< double_<A2>, S2 > ))
                            )
  {
    typedef int result_type;
    
    BOOST_FORCEINLINE result_type operator()(A0& a0, A1 const& a1, A2 const& a2)
    {
      char transa = 'N';
      char transb = 'N';
      long int m = nt2::size(a1)(1); 
      long int n = nt2::size(a2)(2);
      long int k = nt2::size(a1)(2);
      double alpha = 1.0; 
      long int lda = m;
      long int ldb = k;
      double beta = 0.0; 
      long int ldc = m;
      F77NAME(dgemm)(&transa,&transb,&m,&n,&k,&alpha,a1.begin(),&lda,a2.begin(),&ldb,&beta,a0.begin(),&ldc);
    }
  };

/////////////////////////////////////////////////////////////////////////////
// Implementation when type A0 is table_<float_>
/////////////////////////////////////////////////////////////////////////////
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::gemm_, tag::cpu_
                            , (A0)(S0)(A1)(S1)(A2)(S2)
                            , ((table_< single_<A0>, S0 > ))
                              ((table_< single_<A1>, S1 > ))
                              ((table_< single_<A2>, S2 > ))
                            )
  {
    typedef int result_type;

    BOOST_FORCEINLINE result_type operator()(A0& a0, A1 const& a1, A2 const& a2)
    {
      char transa = 'N';
      char transb = 'N'; 
      long int m = nt2::size(a1)(1); 
      long int n = nt2::size(a2)(2);
      long int k = nt2::size(a1)(2);
      float alpha = 1.0; 
      long int lda = m;
      long int ldb = k;
      float beta = 0.0; 
      long int ldc = m;
      F77NAME(sgemm)(&transa,&transb,&m,&n,&k,&alpha,a1.begin(),&lda,a2.begin(),&ldb,&beta,a0.begin(),&ldc);
    }
  };

} }

#endif
