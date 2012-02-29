//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II       
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI       
//                                                                              
//          Distributed under the Boost Software License, Version 1.0.          
//                 See accompanying file LICENSE.txt or copy at                 
//                     http://www.boost.org/LICENSE_1_0.txt                     
//==============================================================================
#ifndef NT2_TOOLBOX_ALGEBRA_FUNCTIONS_LAPACK_LSE_HPP_INCLUDED
#define NT2_TOOLBOX_ALGEBRA_FUNCTIONS_LAPACK_LSE_HPP_INCLUDED

#include <nt2/table.hpp>
#include <nt2/core/settings/forward/shape.hpp>
#include <nt2/core/settings/forward/storage_scheme.hpp>
#include <nt2/core/container/category.hpp>
#include <nt2/include/functions/extent.hpp>
#include <nt2/include/functions/size.hpp>
#include <nt2/toolbox/algebra/lapack/lse.hpp>
#include <nt2/toolbox/algebra/details/padding.hpp>
#include <nt2/include/functions/issquare.hpp>
#include <nt2/include/constants/zero.hpp>
#include <nt2/include/constants/one.hpp>
#include <boost/preprocessor/cat.hpp>
#include <nt2/sdk/error/assert.hpp>
#include <nt2/toolbox/algebra/lapack/lapackworkspace.hpp>
#include <nt2/toolbox/algebra/lapack/l_lse.hpp>
namespace nt2 {

  namespace ext
  {
    //catch all
    NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::lse_, tag::cpu_, 
				(A)(SIZEA)(SHA)(STA)
				(B)(SIZEB)(SHB)(STB)
				(C)(SIZEC)(SHC)(STC)
				(D)(SIZED)(SHD)(STD)
				(X)(SIZEX)(SHX)(STX),  		
				((table_< unspecified_<A>, nt2::settings(SIZEA,SHA,STA) > ))
                                ((table_< unspecified_<B>, nt2::settings(SIZEB,SHB,STB) > ))
                                ((table_< unspecified_<C>, nt2::settings(SIZEC,SHC,STC) > ))
                                ((table_< unspecified_<D>, nt2::settings(SIZED,SHD,STD) > ))
                                ((table_< unspecified_<X>, nt2::settings(SIZEX,SHX,STX) > ))
				)
    {
      typedef void result_type;
      typedef typename C::value_type value_type; 
      
      BOOST_FORCEINLINE result_type operator()(char const&, A&, B&, C&, D& d, X& x)
      {
	BOOST_ASSERT_MSG(false, "The current lapack matrices forms/storages\
                                   is not supported by any lse lapack call");
      }
    };
    
    NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::lse_, tag::cpu_, 
				(A)(SIZEA)
				(B)(SIZEB)
				(C)(SIZEC)
				(D)(SIZED)
				(X)(SIZEX), 
				((table_< floating_<A>, nt2::settings(SIZEA,SHA,STA) > ))
                                ((table_< floating_<B>, nt2::settings(SIZEB,SHB,STB) > ))
                                ((table_< floating_<C>, nt2::settings(SIZEC,SHC,STC) > ))
                                ((table_< floating_<D>, nt2::settings(SIZED,SHD,STD) > ))
                                ((table_< floating_<X>, nt2::settings(SIZEX,SHX,STX) > ))
				
				)
    {
      
      typedef void result_type;
      typedef typename C::value_type value_type; 
      BOOST_FORCEINLINE result_type operator()(A & a, B & b, C& c, D & d, X & x)
      {
	long int m = a.extent(A)[0];
	long int n = a.extent(A)[1];
	long int p = b.extent(A)[0];
	long int lda = nt2::details::padding(a);
	long int ldb = nt2::details::padding(b);
	BOOST_ASSERT_MSG( (n == b.size(2)),
                          "In lse calls the number of columns of a must match the number of rows of b");
	BOOST_ASSERT_MSG( (n >= p),
                          "In lse calls the number of columns of a must be greater than of rows of b");
	BOOST_ASSERT_MSG( (n == x.numel(1)),
                          "In lse calls the number of columns of a must match the number of elements of x)");
	BOOST_ASSERT_MSG( (n <= m+p),
                          "In lse calls n <= m+p");
	BOOST_ASSERT_MSG( (p == d.numel(1)),
                          "In lse calls the number of rows of b must match the number of elements of d");
	//	std::cout << "conventional" << std::endl;
	nt2::details::gelse(&m, &n, &p, a.begin(), &lda, b.begin(), &ldb, c.begin(), d.begin(), &info); 
      }
    };
    
    NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::lse_, tag::cpu_, 
				(A)(SIZEA)
				(B)(SIZEB)
				(C)(SIZEC)
				(D)(SIZED)
				(X)(SIZEX), 
				((table_< complex_<floating_<A> >, nt2::settings(SIZEA,SHA,STA) > ))
                                ((table_< complex_<floating_<B> >, nt2::settings(SIZEB,SHB,STB) > ))
                                ((table_< complex_<floating_<C> >, nt2::settings(SIZEC,SHC,STC) > ))
                                ((table_< complex_<floating_<D> >, nt2::settings(SIZED,SHD,STD) > ))
                                ((table_< complex_<floating_<X> >, nt2::settings(SIZEX,SHX,STX) > ))
				
				)
    {
      
      typedef void result_type;
      typedef typename C::value_type value_type; 
      BOOST_FORCEINLINE result_type operator()(A & a, B & b, C& c, D & d, X & x)
      {
	long int m = a.extent(A)[0];
	long int n = a.extent(A)[1];
	long int p = b.extent(A)[0];
	long int lda = nt2::details::padding(a);
	long int ldb = nt2::details::padding(b);
	BOOST_ASSERT_MSG( (n == b.size(2)),
                          "In lse calls the number of columns of a must match the number of rows of b");
	BOOST_ASSERT_MSG( (n >= p),
                          "In lse calls the number of columns of a must be greater than of rows of b");
	BOOST_ASSERT_MSG( (n == x.numel(1)),
                          "In lse calls the number of columns of a must match the number of elements of x)");
	BOOST_ASSERT_MSG( (n <= m+p),
                          "In lse calls n <= m+p");
	BOOST_ASSERT_MSG( (p == d.numel(1)),
                          "In lse calls the number of rows of b must match the number of elements of d");
	//	std::cout << "conventional" << std::endl;
	nt2::details::gelse(&m, &n, &p, a.begin(), &lda, b.begin(), &ldb, c.begin(), d.begin(), &info);
	//TO DO test info value
      }
    };

  }
}

#endif
