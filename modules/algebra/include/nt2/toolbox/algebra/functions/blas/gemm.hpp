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
#include <nt2/include/functions/size.hpp>
#include <nt2/toolbox/algebra/blas/blas3.hpp>
#include <boost/simd/sdk/memory/align_on.hpp>
#include <nt2/include/constants/zero.hpp>
#include <nt2/include/constants/one.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/assert.hpp>
#include <iostream>

namespace nt2 {
  namespace details {
    template < class T > long int padding(const T & a)
    {
      typedef typename T::parent::lead_t lead_t_a;
      return  boost::simd::memory::align_on(size(a, 1), lead_t_a::value);
    }
  }

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
 
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::gemm_, tag::cpu_, 
                              (A0)(S0)(A1)(S1)(A2)(S2)(A3)(A4)(A5), 
                              ((table_< floating_<A0>, S0 > ))
                              ((table_< floating_<A1>, S1 > ))
                              ((table_< floating_<A2>, S2 > ))
                              (unspecified_ < A5 > ) //TO DO Specify
                              (scalar_ < arithmetic_<A3 > > )
                              (scalar_ < arithmetic_<A4 > > )
                              )
    {
      typedef void result_type;
      typedef A5 targ_t;
      
      BOOST_FORCEINLINE result_type operator()(A0& a0,
                                               A1 const& a1, A2 const& a2,
                                               A5 const& a5, 
                                               A3 const& a3, A4 const& a4
                                               )
      {
        typedef typename A0::value_type value_type; 
        const char transa = targ_t::tA;
        const char transb = targ_t::tB; 
        const long int m = nt2::size(a1)(transa=='T'?2:1); 
        const long int n = nt2::size(a2)(transb=='T'?1:2); 
        const long int k = nt2::size(a1)(transa=='T'?1:2);

        const long int kb = nt2::size(a2)(transb=='T'?2:1);
        BOOST_ASSERT( ( k == kb ) );//inner dimensions must match

        const value_type alpha = a3; 
        const long int lda = nt2::details::padding(a1);
        const long int ldb = nt2::details::padding(a2);
        const value_type beta = a4; 
        const long int ldc = nt2::details::padding(a0);
        gemm(&transa,&transb,&m,&n,&k,&alpha,a1.begin(),&lda,a2.begin(),&ldb,&beta,a0.begin(),&ldc);
      }
    };
  }
  
  template < class T,  class A0,  class A1,  class A2,  class A3,  class A4>
  inline void gemm(A0& a0, A1 const& a1, A2 const& a2,A3 const& a3,A4 const&  a4)
  {
    typedef typename A0::value_type value_type; 
    gemm(a0, a1, a2, T(), a3, a4); 
  }

  template < class T,  class A0,  class A1,  class A2,  class A3>
  inline void gemm(A0& a0, A1 const& a1, A2 const& a2,A3 const& a3)
  {
    typedef typename A0::value_type value_type; 
    gemm(a0, a1, a2, T(), a3, Zero<value_type>()); 
  }
  
  template < class T,  class A0,  class A1,  class A2>
  inline void gemm(A0& a0, A1 const& a1, A2 const& a2)
  {
    typedef typename A0::value_type value_type; 
    gemm(a0, a1, a2, T(), One<value_type>(), Zero<value_type>()); 
  }
}

#endif
