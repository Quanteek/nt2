/*******************************************************************************
 *         Copyright 2003 & onward LASMEA UMR 6602 CNRS/Univ. Clermont II
 *         Copyright 2009 & onward LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_SDK_UNIT_DETAILS_TESTS_HPP_INCLUDED
#define NT2_SDK_UNIT_DETAILS_TESTS_HPP_INCLUDED
#include <nt2/include/functions/ulpdist.hpp>
#include <nt2/sdk/meta/upgrade.hpp>
#include <iostream>

namespace nt2 { namespace details
{
  template<class T, class U, class V>					
  inline bool test_ulp_eq( char const* x1
			   , char const* x2		
			   , char const* x3				
			   , int line
			   , char const * fn			
			   , T const & t
			   , U const & u			
			   , V const & v				
			   )						
  {									
    test_count()++;							
    volatile T tt(t);							
    volatile U uu(u);							
    volatile V vv(v);
    typedef typename nt2::meta::upgrade<T>::type TT;
    typedef typename nt2::meta::upgrade<U>::type UU;
    if( nt2::ulpdist(tt, uu ) <= vv)					
      {									
	std::cout << " * Test `"					
		  << "ulpdist(" << x1 << ", " <<  x2 << ") <= " << x3	
		  << "` **passed**."					
		  << " (" << line << ")"				
		  << std::endl;
	return true; 
      }									
    else								
      {									
	std::cout << " * Test `"					
		  << "ulpdist(" << x1 << ", " <<  x2 << ") <= " << x3	
		  << "` **failed** in function "			
		  << fn << " (" << line << ")"				
		  << "ulpdist(" << TT(t) << ", " <<  UU(u) << ") == "		
		  <<  nt2::ulpdist(tt, uu )				
		  << std::endl;						
	++error_count();
	return false; 
      }									
  }									

} }
#endif
