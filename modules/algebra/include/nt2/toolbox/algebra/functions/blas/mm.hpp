//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II       
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI       
//                                                                              
//          Distributed under the Boost Software License, Version 1.0.          
//                 See accompanying file LICENSE.txt or copy at                 
//                     http://www.boost.org/LICENSE_1_0.txt                     
//==============================================================================
#ifndef NT2_TOOLBOX_ALGEBRA_FUNCTIONS_BLAS_MM_HPP_INCLUDED
#define NT2_TOOLBOX_ALGEBRA_FUNCTIONS_BLAS_MM_HPP_INCLUDED

#include <nt2/table.hpp>
#include <nt2/core/container/category.hpp>
#include <nt2/include/functions/extent.hpp>
#include <nt2/include/functions/size.hpp>
#include <nt2/toolbox/algebra/blas/blas3_b.hpp>
#include <nt2/toolbox/algebra/details/padding.hpp>
#include <nt2/include/functions/issquare.hpp>
#include <nt2/include/constants/zero.hpp>
#include <nt2/include/constants/one.hpp>
#include <boost/preprocessor/cat.hpp>
#include <nt2/sdk/error/assert.hpp>

namespace nt2 {
  namespace ext
  {
    
    NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::mm_, tag::cpu_
                            , (A0)(S0)(A1)(S1)(A2)(S2)(A3)(A4)(STATUS) 
                            , (unspecified_ < STATUS > )
                              ((table_< floating_<A0>, S0 > ))
                              ((table_< floating_<A1>, S1 > ))
                              ((table_< floating_<A2>, S2 > ))
                              (scalar_ < arithmetic_<A3 > > )
                              (scalar_ < arithmetic_<A4 > > )
                
                              )
    {
      typedef void result_type;
      typedef typename A2::value_type value_type; 
      
      BOOST_FORCEINLINE result_type operator()( STATUS const& 
                                              , A0 const& a0, A1 const& a1
                                              , A2& a2
                                              , A3 const& a3, A4 const& a4
                                              )
      {
        blas_call<STATUS::type>::mm_call( STATUS(), a0, a1, a2, a3, a4); 
      }
    private:
//       template <std::size_t btype, class Dummy =  void> struct Blas_call   { static const std::size_t value =  0; };
//       template <class Dummy>                  struct Blas_call<blas_types::general, Dummy>  { static const std::size_t value =  1; };
//       template <class Dummy>                  struct Blas_call<blas_types::symetric, Dummy> { static const std::size_t value =  2; };
//       template <class Dummy>                  struct Blas_call<blas_types::hermitian, Dummy>{ static const std::size_t value =  2; };

      template < std::size_t FORM, class Dummy = void> struct blas_call
      {
        static inline
        void mm_call( STATUS const& 
                      , A0 const& a0, A1 const& a1
                      , A2& a2
                      , A3 const& a3, A4 const& a4)
        {
          BOOST_ASSERT_MSG(false, "The current blas form is not supported by the mm blas call"); 
        }
      };
      
      template <class Dummy> struct blas_call < blas_types::general, Dummy> 
      {
        static inline
        void  mm_call( STATUS const& 
                       , A0 const& a0, A1 const& a1
                       , A2& a2
                       , A3 const& a3, A4 const& a4)
        {
          std::cout << "general" << std::endl; 
          const char transa = STATUS::transa; 
          const char transb = STATUS::transb; 
          const long int m = nt2::extent(a0)[transa=='N'?0:1];
          const long int n = nt2::extent(a1)[transb=='N'?1:0];
          const long int k = nt2::extent(a0)[transa=='N'?1:0];
          
          BOOST_ASSERT_MSG( (k == nt2::size(a1, transb=='N'?1:2)),
                            "In matrix-vector product C = al*A°*B°+ be*C (mm) inner dimensions of A° and B° must match");
          BOOST_ASSERT_MSG( ((m == nt2::size(a2,1))&&(n == nt2::size(a2,2))),
                            "In matrix-vector product C = al*A°*B°+ be*C (mm) outer dimensions of A° and B° must match C ones");
          
          const value_type alpha = a3; 
          const value_type beta = a4; 
          const long int lda = nt2::details::padding(a0);
          const long int ldb = nt2::details::padding(a1);
          const long int ldc = nt2::details::padding(a2);
          nt2::details::gemm(&transa,&transb,&m,&n,&k,&alpha,a0.begin(),&lda,a1.begin(),&ldb,&beta,a2.begin(),&ldc);
        }
      };
      
      template <class Dummy> struct blas_call < blas_types::symetric, Dummy> 
      {
        static inline
        void  mm_call( STATUS const&, 
                       A0 const& a0, A1 const& a1, 
                       A2& a2, 
                       A3 const& a3, A4 const& a4)
        {
          std::cout << "symetric" << std::endl; 
          const char side = STATUS::side; 
          const char uplo = STATUS::uplo; 
          const long int n = nt2::extent(a1)[1];
          const long int m = nt2::extent(a1)[0];
          BOOST_ASSERT_MSG( (is_square(a0)),"In symm calls matrix A must be square");
          
          BOOST_ASSERT_MSG( ( nt2::size(a0,1) == (side=='L')?m:n),
                            "In symm calls inner dimensions, according to side, must match");
          
          const value_type alpha = a3; 
          const value_type beta  = a4; 
          const long int lda = nt2::details::padding(a0);
          const long int ldb = nt2::details::padding(a1);
          const long int ldc = nt2::details::padding(a2);
          nt2::details::symm(&side,&uplo,&m,&n,&alpha,a0.begin(),&lda,a1.begin(),&ldb,&beta,a2.begin(),&ldc);
        }
      };
      
      template <class Dummy> struct blas_call < blas_types::hermitian, Dummy> 
      {
        static inline
        void  mm_call( STATUS const&, 
                       A0 const& a0, A1 const& a1, 
                       A2& a2, 
                       A3 const& a3, A4 const& a4)
        {
          std::cout << "hermitian" << std::endl; 
          const char side = STATUS::side; 
          const char uplo = STATUS::uplo; 
          const long int n = nt2::extent(a1)[1];
          const long int m = nt2::extent(a1)[0];
          BOOST_ASSERT_MSG( (is_square(a0)),"In hemm calls matrix A must be square");
          
          BOOST_ASSERT_MSG( ( nt2::size(a0,1) == (side=='L')?m:n) ,
                            "In hemm calls, inner dimensions according to side, must match");
          const value_type alpha = a3; 
          const long int lda = nt2::details::padding(a0);
          const long int ldb = nt2::details::padding(a1);
          const value_type beta = a4; 
          const long int ldc = nt2::details::padding(a2);
          nt2::details::hemm(&side,&uplo,&m,&n,&alpha,a0.begin(),&lda,a1.begin(),&ldb,&beta,a2.begin(),&ldc);
        }
      }; 
    };
  }
}

#endif
