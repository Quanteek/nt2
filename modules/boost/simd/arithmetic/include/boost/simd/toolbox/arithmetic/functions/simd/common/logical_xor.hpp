//==============================================================================
//         Copyright 2003 - 2011 LASMEA UMR 6602 CNRS/Univ. Clermont II         
//         Copyright 2009 - 2011 LRI    UMR 8623 CNRS/Univ Paris Sud XI         
//                                                                              
//          Distributed under the Boost Software License, Version 1.0.          
//                 See accompanying file LICENSE.txt or copy at                 
//                     http://www.boost.org/LICENSE_1_0.txt                     
//==============================================================================
#ifndef BOOST_SIMD_TOOLBOX_ARITHMETIC_FUNCTIONS_SIMD_COMMON_LOGICAL_XOR_HPP_INCLUDED
#define BOOST_SIMD_TOOLBOX_ARITHMETIC_FUNCTIONS_SIMD_COMMON_LOGICAL_XOR_HPP_INCLUDED
#include <boost/simd/include/functions/is_nez.hpp>
#include <boost/simd/include/functions/bitwise_xor.hpp>
#include <boost/simd/include/functions/genmask.hpp>
// TODO put it in boolean or operator ?

namespace boost { namespace simd { namespace ext
{
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( boost::simd::tag::logical_xor_, tag::cpu_, (A0)(X)
                            , ((simd_<arithmetic_<A0>,X>))
                              ((simd_<arithmetic_<A0>,X>))
                            )
  {
    typedef typename meta::as_logical<A0>::type result_type;
    BOOST_SIMD_FUNCTOR_CALL_REPEAT(2)
      {
       return is_nez(b_xor(genmask<A0>(a0),genmask<A0>(a1)));
      }
  };
  
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( boost::simd::tag::logical_xor_, tag::cpu_, (A0)(X)
                            , ((simd_<logical_<A0>,X>))
                              ((simd_<logical_<A0>,X>))
                            )
  {
    typedef typename meta::as_logical<A0>::type result_type;
    BOOST_SIMD_FUNCTOR_CALL_REPEAT(2)
      {
       typedef typename A0::type type; 
       return is_nez(b_xor(genmask(a0),genmask(a1)));
      }
  };
} } }

#endif
