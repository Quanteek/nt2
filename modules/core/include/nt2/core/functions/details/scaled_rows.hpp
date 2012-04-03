//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#ifndef NT2_CORE_FUNCTIONS_DETAILS_SCALED_ROWS_HPP_INCLUDED
#define NT2_CORE_FUNCTIONS_DETAILS_SCALED_ROWS_HPP_INCLUDED

#include <boost/fusion/include/at.hpp>
#include <nt2/include/functions/if_else.hpp>
#include <nt2/include/functions/is_equal.hpp>
#include <nt2/include/functions/splat.hpp>
#include <nt2/include/functions/enumerate.hpp>
#include <nt2/include/functions/arith.hpp>
#include <nt2/include/constants/one.hpp>
#include <iostream>

namespace nt2 { namespace details
{

  template < class T, class T1>
  struct scaled_rows
  {
    scaled_rows()                : start_(T()), h_(One<T>())  {}
    scaled_rows(const T & start, const T1 & h) : start_(start), h_(h){}
    template<class Pos, class Size, class Target>
    typename Target::type operator()(Pos const& p, Size const&, Target const& ) const
    {
      const  int b = 1; 
      static int i =  0; 
      typedef typename Target::type result_type;
      typedef typename meta::scalar_of<result_type>::type s_type;
      std::cout << "sizeof " << sizeof(result_type) << " i    -> " << ++i << "p " << boost::fusion::at_c<0>(p) << std::endl;
      std::cout << "result (" << i << "--   "<<  nt2::arith<result_type>(boost::fusion::at_c<0>(p)+s_type(start_-b), s_type(h_)) <<"  --" << i<< ")" << std::endl;
      return nt2::arith<result_type>(boost::fusion::at_c<0>(p)+s_type(start_-b), s_type(h_));
    }
  private :
    T start_;
    T h_; 
  };
} }

#endif
