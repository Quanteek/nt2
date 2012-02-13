//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II       
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI       
//                                                                              
//          Distributed under the Boost Software License, Version 1.0.          
//                 See accompanying file LICENSE.txt or copy at                 
//                     http://www.boost.org/LICENSE_1_0.txt                     
//==============================================================================
#ifndef NT2_TOOLBOX_ALGEBRA_FUNCTIONS_BLAS_GENERAL_GEMM_HPP_INCLUDED
#define NT2_TOOLBOX_ALGEBRA_FUNCTIONS_BLAS_GENERAL_GEMM_HPP_INCLUDED

#include <nt2/table.hpp>
#include <nt2/core/container/category.hpp>
#include <nt2/include/functions/size.hpp>
#include <nt2/toolbox/algebra/blas/blas3.hpp>
#include <boost/simd/sdk/memory/align_on.hpp>

#define F77NAME(x) x##_

namespace nt2 {
  struct NN { static const char tA = 'N'; static const char tB = 'N'; };
  struct TN { static const char tA = 'T'; static const char tB = 'N'; };
  struct NT { static const char tA = 'N'; static const char tB = 'T'; };
  struct TT { static const char tA = 'T'; static const char tB = 'T'; };
  
  namespace ext
  {
    
#define NT2_GEMM(T, PREFIX)                                             \
    inline void gemm(const char *ta, const char *tb, const long int *m, \
                     const long int *n, const long int *k,              \
                     const T *al, const T *a,                           \
                     const long int *lda, const T *b,                   \
                     const long int *ldb, const T *be, T *c,            \
                     const long int *ldc)                               \
    {                                                                   \
      F77NAME(PREFIX ## gemm)(ta,tb,m,n,k,al,a,lda,b,ldb,be,c,ldc);     \
    }                                                                   \

    NT2_GEMM(double, d)
    NT2_GEMM(float,  s)
    NT2_GEMM(std::complex<double>, z)
    NT2_GEMM(std::complex<float>, c)

#undef NT2_GEMM
 
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::general_gemm_, tag::cpu_, 
                              (A0)(S0)(A1)(S1)(A2)(S2)(A3)(A4)(A5), 
                              ((table_< floating_<A0>, S0 > ))
                              ((table_< floating_<A1>, S1 > ))
                              ((table_< floating_<A2>, S2 > ))
                              (unspecified_ < A5 > )
                              (scalar_ < arithmetic_<A3 > > )
                              (scalar_ < arithmetic_<A4 > > )

                            )
  {
    typedef int result_type;
    typedef typename A0::parent::lead_t lead_t_a0;
    typedef typename A1::parent::lead_t lead_t_a1;
    typedef typename A2::parent::lead_t lead_t_a2;
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
      const long int n = nt2::size(a2)(transa=='T'?1:2); 
      const long int k = nt2::size(a1)(transa=='T'?1:2);
      std::cout << "m =  " << m
                << " n =  " << n
                << " k =  " << k << std::endl; 
      const value_type alpha = a3; 
      const long int lda = boost::simd::memory::align_on(m, lead_t_a1::value);
      const long int ldb = boost::simd::memory::align_on(k, lead_t_a2::value);
      const value_type beta = a4; 
      const long int ldc = boost::simd::memory::align_on(m, lead_t_a1::value);
      std::cout << "lda =  " << lda
                << " ldb =  " << ldb
                << " ldc =  " << ldc << std::endl; 
      gemm(&transa,&transb,&m,&n,&k,&alpha,a1.begin(),&lda,a2.begin(),&ldb,&beta,a0.begin(),&ldc);
      return 0; 
    }
   };

// /////////////////////////////////////////////////////////////////////////////
// // Implementation when type A0 is table_<float_>
// /////////////////////////////////////////////////////////////////////////////
//   NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::gemm_, tag::cpu_
//                             , (A0)(S0)(A1)(S1)(A2)(S2)
//                             , ((table_< single_<A0>, S0 > ))
//                               ((table_< single_<A1>, S1 > ))
//                               ((table_< single_<A2>, S2 > ))
//                               (scalar_ < A3 > )
//                               (scalar_ < A4 > )
//                             )
//   {
//     typedef void result_type;
//     typedef typename A1::parent::lead_t lead_t_a1;
//     typedef typename A2::parent::lead_t lead_t_a2;
    
//     BOOST_FORCEINLINE result_type operator()(A0& a0, A1 const& a1, A2 const& a2,
//                                                      A3 const& a3 = 1, A4 const& a4 = 0)
//     {
//       const char transa = 'N';
//       const char transb = 'N';
//       const long int m = nt2::size(a1)(1);
//       const long int n = nt2::size(a2)(2);
//       const long int k = nt2::size(a1)(2);
//       const float alpha = a3; 
//       const long int lda = boost::simd::memory::align_on(m, lead_t_a1::value);
//       const long int ldb = boost::simd::memory::align_on(k, lead_t_a2::value);
//       const float beta = a4; 
//       const long int ldc = boost::simd::memory::align_on(m, lead_t_a1::value);
//       gemm(&transa,&transb,&m,&n,&k,&alpha,a1.begin(),&lda,a2.begin(),&ldb,&beta,a0.begin(),&ldc);
//     }
   
//   };

} }

#undef F77NAME

#endif
