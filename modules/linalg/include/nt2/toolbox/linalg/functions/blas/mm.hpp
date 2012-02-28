//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II       
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI       
//                                                                              
//          Distributed under the Boost Software License, Version 1.0.          
//                 See accompanying file LICENSE.txt or copy at                 
//                     http://www.boost.org/LICENSE_1_0.txt                     
//==============================================================================
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_BLAS_MM_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_BLAS_MM_HPP_INCLUDED

#include <nt2/table.hpp>
#include <nt2/core/settings/forward/shape.hpp>
#include <nt2/core/settings/forward/storage_scheme.hpp>
#include <nt2/core/container/category.hpp>
#include <nt2/include/functions/extent.hpp>
#include <nt2/include/functions/size.hpp>
#include <nt2/toolbox/linalg/details/blas/mm.hpp>
#include <nt2/toolbox/linalg/details/utility/padding.hpp>
#include <nt2/include/functions/issquare.hpp>
#include <nt2/include/constants/zero.hpp>
#include <nt2/include/constants/one.hpp>
#include <boost/preprocessor/cat.hpp>
#include <nt2/sdk/error/assert.hpp>

namespace nt2 {
  template < class T > T trans(T const & a){ return a; }// this is fake transpose
  template < class T > T transconj(T const & a){ return a; }// this is fake tranconjugate
  namespace details
  {
    template < class T > struct run_status
    {
      long int m, n, k;
      long int lda, ldb, ldc;
      T alpha, beta; 
      char uplo;
      char transa;
      char transb; 
      char diag;
      char side;
    };
    
    template < class STATUS, class A, class B, class C, class Alpha, class Beta, class RSTATUS> 
    inline void mm_prepare(STATUS const&, A const& a, B const& b, C const& c,
                           Alpha const& alpha, Beta const& beta, RSTATUS & r)
      {
        r.transa = STATUS::transa; 
        r.transb = STATUS::transb;
        //        r.side  =  STATUS::side;
        //        r.diag  =  STATUS::diag;
        r.uplo  =  STATUS::uplo; 
        r.m = nt2::extent(a)[r.transa=='N'?0:1];
        r.n = nt2::extent(b)[r.transb=='N'?1:0];
        r.k = nt2::extent(a)[r.transa=='N'?1:0];
        r.alpha = alpha; 
        r.beta = beta; 
        r.lda = nt2::details::padding(a);
        r.ldb = nt2::details::padding(b);
        r.ldc = nt2::details::padding(c);
      }
  }
  namespace ext
  {
    
    NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::mm_, tag::cpu_, 
                             (A)(SIZEA)(SHA)(STA)
                             (B)(SIZEB)(SHB)(STB)
                             (C)(SIZEC)(SHC)(STC)(Alpha)(Beta)(STATUS), 
                             (unspecified_ < STATUS > )
                                ((table_< floating_<A>, nt2::settings(SIZEA,SHA,STA) > ))
                                ((table_< floating_<B>, nt2::settings(SIZEB,SHB,STB) > ))
                                ((table_< floating_<C>, nt2::settings(SIZEC,SHC,STC) > ))
                             (scalar_ < arithmetic_<Alpha > > )
                             (scalar_ < arithmetic_<Beta > > )
                              )
    {
      typedef void result_type;
      typedef typename C::value_type value_type; 
      
      BOOST_FORCEINLINE result_type operator()( STATUS const& 
                                              , A const& a, B const& b
                                              , C& c
                                              , Alpha const& alpha, Beta const& beta
                                              )
      {
          BOOST_ASSERT_MSG(false, "The current blas matrices forms/storages\
                                   is not supported by any mm blas call");
      }
    };

    NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::mm_, tag::cpu_, 
                             (A)(SIZEA)
                             (B)(SIZEB)
                             (C)(SIZEC)(Alpha)(Beta)(STATUS), 
                             (unspecified_ < STATUS > )
                                ((table_< floating_<A>, nt2::settings(SIZEA,nt2::rectangular_,nt2::conventional_) > ))
                                ((table_< floating_<B>, nt2::settings(SIZEB,nt2::rectangular_,nt2::conventional_) > ))
                                ((table_< floating_<C>, nt2::settings(SIZEC,nt2::rectangular_,nt2::conventional_) > ))
                             (scalar_ < arithmetic_<Alpha > > )
                             (scalar_ < arithmetic_<Beta > > )
                              )
    {
      typedef void result_type;
      typedef typename C::value_type value_type; 
      BOOST_FORCEINLINE result_type operator()( STATUS const&, 
                                                A const& a, B const& b, 
                                                C& c, 
                                                Alpha const& alpha, Beta const& beta
                                                )
      {
        std::cout << "conventional" << std::endl;
        details::run_status<value_type> r; 
        details::mm_prepare(STATUS(),a,b,c,alpha,beta,r);
        BOOST_ASSERT_MSG( (r.k == nt2::size(b, r.transb=='N'?1:2)),
                          "In matrix-vector product C = al*A°*B°+ be*C (gemm) inner dimensions of A° and B° must match");
        BOOST_ASSERT_MSG( ((r.m == nt2::size(c,1))&&(r.n == nt2::size(c,2))),
                          "In matrix-vector product C = al*A°*B°+ be*C (gemm) outer dimensions of A° and B° must match C ones");
        nt2::details::gemm(&r.transa,&r.transb,&r.m,&r.n,&r.k,&r.alpha,a.begin(),
                           &r.lda,b.begin(),&r.ldb,&r.beta,c.begin(),&r.ldc);
      }
    };

    NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::mm_, tag::cpu_, 
                             (A)(SIZEA)
                             (B)(SIZEB)
                             (C)(SIZEC)(Alpha)(Beta)(STATUS), 
                             (unspecified_ < STATUS > )
                                ((table_< complex_<floating_<A> >, nt2::settings(SIZEA,nt2::rectangular_,nt2::conventional_) > ))
                                ((table_< complex_<floating_<B> >, nt2::settings(SIZEB,nt2::rectangular_,nt2::conventional_) > ))
                                ((table_< complex_<floating_<C> >, nt2::settings(SIZEC,nt2::rectangular_,nt2::conventional_) > ))
                             (scalar_ < complex_<arithmetic_<Alpha > > > )
                             (scalar_ < complex_<arithmetic_<Beta > > > )
                              )
    {
      typedef void result_type;
      typedef typename C::value_type value_type; 
      BOOST_FORCEINLINE result_type operator()( STATUS const&, 
                                                A const& a, B const& b, 
                                                C& c, 
                                                Alpha const& alpha, Beta const& beta
                                                )
      {
        std::cout << "conventional" << std::endl;
        details::run_status<value_type> r; 
        details::mm_prepare(STATUS(),a,b,c,alpha,beta,r);
        BOOST_ASSERT_MSG( (r.k == nt2::size(b, r.transb=='N'?1:2)),
                          "In matrix-vector product C = al*A°*B°+ be*C (gemm) inner dimensions of A° and B° must match");
        BOOST_ASSERT_MSG( ((r.m == nt2::size(c,1))&&(r.n == nt2::size(c,2))),
                          "In matrix-vector product C = al*A°*B°+ be*C (gemm) outer dimensions of A° and B° must match C ones");
        nt2::details::gemm(&r.transa,&r.transb,&r.m,&r.n,&r.k,&r.alpha,a.begin(),
                           &r.lda,b.begin(),&r.ldb,&r.beta,c.begin(),&r.ldc);
      }
    };
    
    NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::mm_, tag::cpu_, 
                                (A)(SIZEA)
                                (B)(SIZEB)
                                (C)(SIZEC)(Alpha)(Beta)(STATUS), 
                                (unspecified_ < STATUS > )
                                ((table_< floating_<A>, nt2::settings(SIZEA,nt2::symetric_,nt2::conventional_) > ))
                                ((table_< floating_<B>, nt2::settings(SIZEB,nt2::rectangular_,nt2::conventional_) > ))
                                ((table_< floating_<C>, nt2::settings(SIZEC,nt2::rectangular_,nt2::conventional_) > ))
                                (scalar_ < arithmetic_<Alpha > > )
                                (scalar_ < arithmetic_<Beta > > )
                                )
    {
      typedef void result_type;
      typedef typename C::value_type value_type; 
      
      BOOST_FORCEINLINE result_type operator()( STATUS const& 
                                                , A const& a, B const& b
                                                , C& c
                                                , Alpha const& alpha, Beta const& beta
                                                )
      {
        std::cout << "symetric a" << std::endl;
        details::run_status<value_type> r; 
        details::mm_prepare(STATUS(),a,b,c,alpha,beta,r);
//         BOOST_ASSERT_MSG( (is_square(a)),"In symm calls matrix A must be square");//Perhaps no use as we know A is symetric_
//         BOOST_ASSERT_MSG( ( nt2::size(a,1) == r.m)  "In mm calls inner dimensions, must match"); 
//         BOOST_ASSERT_MSG( ( nt2::size(a,1) == (r.side=='L')?r.m:r.n),
//                           "In symm calls inner dimensions, according to side, must match");
        if (r.transb != 'N')
          {
            const long int side = 'L'; 
            nt2::details::symm(side,&r.uplo,&r.m,&r.n,&r.alpha,a.begin(),
                               &r.lda,b.begin(),&r.ldb,&r.beta,c.begin(),&r.ldc);
          }
        else
          {
            const long int side = 'R'; 
            C c1 = nt2::trans(c);  
            const long int ldc1 = nt2::details::padding(c1);
            nt2::details::symm(side,&r.uplo,&r.m,&r.n,&r.alpha,a.begin(),
                               &r.lda,b.begin(),&r.ldb,&r.beta,c1.begin(),ldc1);
            c =  nt2::trans(c1); 
          }
        
      }
    };

    NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::mm_, tag::cpu_, 
                                (A)(SIZEA)
                                (B)(SIZEB)
                                (C)(SIZEC)(Alpha)(Beta)(STATUS), 
                                (unspecified_ < STATUS > )
                                ((table_< floating_<A>, nt2::settings(SIZEA,nt2::rectangular_,nt2::conventional_) > ))
                                ((table_< floating_<B>, nt2::settings(SIZEB,nt2::symetric_,nt2::conventional_) > ))
                                ((table_< floating_<C>, nt2::settings(SIZEC,nt2::rectangular_,nt2::conventional_) > ))
                                (scalar_ < arithmetic_<Alpha > > )
                                (scalar_ < arithmetic_<Beta > > )
                                )
    {
      typedef void result_type;
      typedef typename C::value_type value_type; 
      
      BOOST_FORCEINLINE result_type operator()( STATUS const& 
                                                , A const& a, B const& b
                                                , C& c
                                                , Alpha const& alpha, Beta const& beta
                                                )
      {
        std::cout << "symetric b" << std::endl;
        details::run_status<value_type> r; 
        details::mm_prepare(STATUS(),a,b,c,alpha,beta,r);
//         BOOST_ASSERT_MSG( (is_square(b)),"In symm calls matrix A must be square");//Perhaps no use as we know A is symetric_
//         BOOST_ASSERT_MSG( ( nt2::size(a,1) == r.m)); 
        //         BOOST_ASSERT_MSG( ( nt2::size(a,1) == (r.side=='L')?r.m:r.n),
        //                           "In symm calls inner dimensions, according to side, must match");
        if (r.transb != 'N')
          {
            const long int side = 'R'; 
            nt2::details::symm(side,&r.uplo,&r.m,&r.n,&r.alpha,a.begin(),
                               &r.lda,b.begin(),&r.ldb,&r.beta,c.begin(),&r.ldc);
          }
        else
          {
            const long int side = 'L'; 
            C c1 = nt2::trans(c);  
            const long int ldc1 = nt2::details::padding(c1);
            nt2::details::symm(side,&r.uplo,&r.m,&r.n,&r.alpha,a.begin(),
                               &r.lda,b.begin(),&r.ldb,&r.beta,c1.begin(),ldc1);
            c =  nt2::trans(c1); 
          }
        
      }
    };
    NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::mm_, tag::cpu_, 
                                (A)(SIZEA)
                                (B)(SIZEB)
                                (C)(SIZEC)(Alpha)(Beta)(STATUS), 
                                (unspecified_ < STATUS > )
                                ((table_< complex_<floating_<A> >, nt2::settings(SIZEA,nt2::symetric_,nt2::conventional_) > ))
                                ((table_< complex_<floating_<B> >, nt2::settings(SIZEB,nt2::rectangular_,nt2::conventional_) > ))
                                ((table_< complex_<floating_<C> >, nt2::settings(SIZEC,nt2::rectangular_,nt2::conventional_) > ))
                                (scalar_ < complex_<arithmetic_<Alpha > > > )
                                (scalar_ < complex_<arithmetic_<Beta > > > )
                                )
    {
      typedef void result_type;
      typedef typename C::value_type value_type; 
      
      BOOST_FORCEINLINE result_type operator()( STATUS const& 
                                                , A const& a, B const& b
                                                , C& c
                                                , Alpha const& alpha, Beta const& beta
                                                )
      {
        std::cout << "symetric a" << std::endl;
        details::run_status<value_type> r; 
        details::mm_prepare(STATUS(),a,b,c,alpha,beta,r);
//         BOOST_ASSERT_MSG( (is_square(a)),"In symm calls matrix A must be square");//Perhaps no use as we know A is symetric_
//         BOOST_ASSERT_MSG( ( nt2::size(a,1) == r.m)  "In mm calls inner dimensions, must match"); 
//         BOOST_ASSERT_MSG( ( nt2::size(a,1) == (r.side=='L')?r.m:r.n),
//                           "In symm calls inner dimensions, according to side, must match");
        if (r.transb != 'N')
          {
            const long int side = 'L'; 
            nt2::details::symm(side,&r.uplo,&r.m,&r.n,&r.alpha,a.begin(),
                               &r.lda,b.begin(),&r.ldb,&r.beta,c.begin(),&r.ldc);
          }
        else
          {
            const long int side = 'R'; 
            C c1 = nt2::trans(c);  
            const long int ldc1 = nt2::details::padding(c1);
            nt2::details::symm(side,&r.uplo,&r.m,&r.n,&r.alpha,a.begin(),
                               &r.lda,b.begin(),&r.ldb,&r.beta,c1.begin(),ldc1);
            c =  nt2::trans(c1); 
          }
        
      }
    };

    NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::mm_, tag::cpu_, 
                                (A)(SIZEA)
                                (B)(SIZEB)
                                (C)(SIZEC)(Alpha)(Beta)(STATUS), 
                                (unspecified_ < STATUS > )
                                ((table_< complex_<floating_<A> >, nt2::settings(SIZEA,nt2::rectangular_,nt2::conventional_) > ))
                                ((table_< complex_<floating_<B> >, nt2::settings(SIZEB,nt2::symetric_,nt2::conventional_) > ))
                                ((table_< complex_<floating_<C> >, nt2::settings(SIZEC,nt2::rectangular_,nt2::conventional_) > ))
                                (scalar_ < complex_<arithmetic_<Alpha > > > )
                                (scalar_ < complex_<arithmetic_<Beta > > > )
                                )
    {
      typedef void result_type;
      typedef typename C::value_type value_type; 
      
      BOOST_FORCEINLINE result_type operator()( STATUS const& 
                                                , A const& a, B const& b
                                                , C& c
                                                , Alpha const& alpha, Beta const& beta
                                                )
      {
        std::cout << "symetric b" << std::endl;
        details::run_status<value_type> r; 
        details::mm_prepare(STATUS(),a,b,c,alpha,beta,r);
//         BOOST_ASSERT_MSG( (is_square(b)),"In symm calls matrix A must be square");//Perhaps no use as we know A is symetric_
//         BOOST_ASSERT_MSG( ( nt2::size(a,1) == r.m)); 
        //         BOOST_ASSERT_MSG( ( nt2::size(a,1) == (r.side=='L')?r.m:r.n),
        //                           "In symm calls inner dimensions, according to side, must match");
        if (r.transb != 'N')
          {
            const long int side = 'R'; 
            nt2::details::symm(side,&r.uplo,&r.m,&r.n,&r.alpha,a.begin(),
                               &r.lda,b.begin(),&r.ldb,&r.beta,c.begin(),&r.ldc);
          }
        else
          {
            const long int side = 'L'; 
            C c1 = nt2::trans(c);  
            const long int ldc1 = nt2::details::padding(c1);
            nt2::details::symm(side,&r.uplo,&r.m,&r.n,&r.alpha,a.begin(),
                               &r.lda,b.begin(),&r.ldb,&r.beta,c1.begin(),ldc1);
            c =  nt2::trans(c1); 
          }
        
      }
    };

    NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::mm_, tag::cpu_, 
                                (A)(SIZEA)
                                (B)(SIZEB)
                                (C)(SIZEC)(Alpha)(Beta)(STATUS), 
                                (unspecified_ < STATUS > )
                                ((table_< complex_<floating_<A> >, nt2::settings(SIZEA,nt2::hermitian_,nt2::conventional_) > ))
                                ((table_< complex_<floating_<B> >, nt2::settings(SIZEB,nt2::rectangular_,nt2::conventional_) > ))
                                ((table_< complex_<floating_<C> >, nt2::settings(SIZEC,nt2::rectangular_,nt2::conventional_) > ))
                                (scalar_< complex_< arithmetic_<Alpha > > > )
                                (scalar_< complex_< arithmetic_<Beta > > > )
                                )
    {
      typedef void result_type;
      typedef typename C::value_type value_type; 
      
      BOOST_FORCEINLINE result_type operator()( STATUS const& 
                                                , A const& a, B const& b
                                                , C& c
                                                , Alpha const& alpha, Beta const& beta
                                                )
      {
        std::cout << "hermitian a" << std::endl;
        details::run_status<value_type> r; 
        details::mm_prepare(STATUS(),a,b,c,alpha,beta,r);
//         BOOST_ASSERT_MSG( (is_square(a)),"In hemm calls matrix A must be square");//Perhaps no use as we know A is hermitian_
//         BOOST_ASSERT_MSG( ( nt2::size(a,1) == r.m)); 
//         BOOST_ASSERT_MSG( ( nt2::size(a,1) == (r.side=='L')?r.m:r.n),
//                           "In hemm calls inner dimensions, according to side, must match");
        if (r.transb != 'N')
          {
            const long int side = 'L'; 
            nt2::details::hemm(side,&r.uplo,&r.m,&r.n,&r.alpha,a.begin(),
                               &r.lda,b.begin(),&r.ldb,&r.beta,c.begin(),&r.ldc);
          }
        else
          {
            const long int side = 'R'; 
            C c1 = nt2::transconj(c);  
            const long int ldc1 = nt2::details::padding(c1);
            nt2::details::hemm(side,&r.uplo,&r.m,&r.n,&r.alpha,a.begin(),
                               &r.lda,b.begin(),&r.ldb,&r.beta,c1.begin(),ldc1);
            c =  nt2::transconj(c1); 
          }
        BOOST_ASSERT_MSG( (is_square(a)),"In hemm calls matrix A must be square");
        
      }
    };

    NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::mm_, tag::cpu_, 
                                (A)(SIZEA)
                                (B)(SIZEB)
                                (C)(SIZEC)(Alpha)(Beta)(STATUS), 
                                (unspecified_ < STATUS > )
                                ((table_< complex_<floating_<A> >, nt2::settings(SIZEA,nt2::rectangular_,nt2::conventional_) > ))
                                ((table_< complex_<floating_<B> >, nt2::settings(SIZEB,nt2::hermitian_,nt2::conventional_) > ))
                                ((table_< complex_<floating_<C> >, nt2::settings(SIZEC,nt2::rectangular_,nt2::conventional_) > ))
                                (scalar_< complex_< arithmetic_<Alpha > > > )
                                (scalar_< complex_< arithmetic_<Beta > > > )
                                )
    {
      typedef void result_type;
      typedef typename C::value_type value_type; 
      
      BOOST_FORCEINLINE result_type operator()( STATUS const& 
                                                , A const& a, B const& b
                                                , C& c
                                                , Alpha const& alpha, Beta const& beta
                                                )
      {
        std::cout << "hermitian b" << std::endl;
        details::run_status<value_type> r; 
        details::mm_prepare(STATUS(),a,b,c,alpha,beta,r);
//         BOOST_ASSERT_MSG( (is_square(b)),"In hemm calls matrix A must be square");//Perhaps no use as we know A is hermitian_
//         BOOST_ASSERT_MSG( ( nt2::size(a,1) == r.m)); 
        //         BOOST_ASSERT_MSG( ( nt2::size(a,1) == (r.side=='L')?r.m:r.n),
        //                           "In hemm calls inner dimensions, according to side, must match");
        if (r.transb != 'N')
          {
            const long int side = 'R'; 
            nt2::details::hemm(side,&r.uplo,&r.m,&r.n,&r.alpha,a.begin(),
                               &r.lda,b.begin(),&r.ldb,&r.beta,c.begin(),&r.ldc);
          }
        else
          {
            const long int side = 'L'; 
            C c1 = nt2::transconj(c);  
            const long int ldc1 = nt2::details::padding(c1);
            nt2::details::hemm(side,&r.uplo,&r.m,&r.n,&r.alpha,a.begin(),
                               &r.lda,b.begin(),&r.ldb,&r.beta,c1.begin(),ldc1);
            c =  nt2::transconj(c1); 
          }
        
      }
    };

  }
}

#endif
