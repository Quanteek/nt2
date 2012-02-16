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
#include <nt2/core/container/category.hpp>
#include <nt2/include/functions/extent.hpp>
#include <nt2/include/functions/size.hpp>
#include <nt2/toolbox/algebra/blas/blas3.hpp>
#include <nt2/toolbox/algebra/details/padding.hpp>
#include <nt2/include/constants/zero.hpp>
#include <nt2/include/constants/one.hpp>
#include <boost/preprocessor/cat.hpp>
#include <nt2/sdk/error/assert.hpp>

namespace nt2 {

  template<char T0, char T1>
  struct gemm_status
  {
    static const char tA = T0;
    static const char tB = T1;
  };

  namespace ext
  {
    
    #define NT2_GEMM(T, PREFIX)                                         \
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

    NT2_GEMM(double, d)
    NT2_GEMM(float,  s)
    NT2_GEMM(std::complex<double>, z)
    NT2_GEMM(std::complex<float>, c)

#undef NT2_GEMM
 
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::gemm_, tag::cpu_
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
                                              , A0 const& a0, A1 const& a1
                                              , A2& a2
                                              , A3 const& a3, A4 const& a4
                                              )
      {
        typedef typename A0::value_type value_type; 

        const char transa = A5::tA; 
        const char transb = A5::tB; 
        const long int m = nt2::extent(a0)[transa=='N'?0:1]; // nt2::size(a0)(transa=='T'?2:1); 
        const long int n = nt2::extent(a1)[transb=='N'?1:0]; // nt2::size(a1)(transb=='T'?1:2); 
        const long int k = nt2::extent(a0)[transa=='N'?1:0]; // nt2::size(a0)(transa=='T'?1:2);

        BOOST_ASSERT_MSG( (k == nt2::size(a1, transb=='N'?1:2)),
                          "In matrix-vector product C = al*A*B+ be*C (gemm) inner dimensions of A and B must match");
                          
                          const value_type alpha = a3; 
        const long int lda = nt2::details::padding(a0);
        const long int ldb = nt2::details::padding(a1);
        const value_type beta = a4; 
        const long int ldc = nt2::details::padding(a2);
        gemm(&transa,&transb,&m,&n,&k,&alpha,a0.begin(),&lda,a1.begin(),&ldb,&beta,a2.begin(),&ldc);
      }
    };
  }
  //call with status
  template < class A5,  class A0,  class A1,  class A2,  class A3>
  inline void gemm(A5 const& a5, A0 const& a0, A1 const& a1, A2& a2, A3 const& a3)
  {
    typedef typename A0::value_type value_type; 
    gemm(a5, a0, a1, a2, a3, Zero<value_type>()); 
  }
  
  template < class A5,  class A0,  class A1,  class A2>
  inline void gemm(A5 const& a5, A0 const& a0, A1 const& a1, A2& a2)
  {
    typedef typename A0::value_type value_type; 
    gemm(a5, a0, a1, a2, One<value_type>(), Zero<value_type>()); 
  }

  //templated on status
  template < char transa,  char transb,  class A0,  class A1,  class A2,  class A3,  class A4>
  inline void gemm(A0 const& a0, A1 const& a1, A2& a2, A3 const& a3, A4 const& a4)
  {
    typedef typename A0::value_type value_type; 
    gemm(gemm_status<transa, transb>(), a0, a1, a2, a3, a4); 
  }

  template < char transa,  char transb,  class A0,  class A1,  class A2,  class A3>
  inline void gemm(A0 const& a0, A1 const& a1, A2& a2, A3 const& a3)
  {
    typedef typename A0::value_type value_type; 
    gemm(gemm_status<transa, transb>(), a0, a1, a2, a3, Zero<value_type>()); 
  }
  
  template < char transa,  char transb, class A0,  class A1,  class A2>
  inline void gemm(A0 const& a0, A1 const& a1, A2& a2)
  {
    typedef typename A0::value_type value_type; 
    gemm(gemm_status<transa, transb>(), a0, a1, a2, One<value_type>(), Zero<value_type>()); 
  }
}

#endif
