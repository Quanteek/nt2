//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#ifndef BOOST_SIMD_TOOLBOX_CONSTANT_CONSTANTS_DETAILS_DIGITS_HPP_INCLUDED
#define BOOST_SIMD_TOOLBOX_CONSTANT_CONSTANTS_DETAILS_DIGITS_HPP_INCLUDED

#include <boost/simd/include/functions/splat.hpp>

//==============================================================================
// Register dispatch over digits<N>
//==============================================================================
namespace boost { namespace dispatch { namespace meta
{
  BOOST_DISPATCH_FUNCTOR_IMPLEMENTATION_TPL( boost::simd::tag::digit_<N> , tag::cpu_
                                , (boost::simd::int64_t N)(class A0)
                                , (target_< scalar_< fundamental_<A0> > >)
                                )
  {
    typedef typename A0::type result_type;

    BOOST_DISPATCH_FUNCTOR_CALL(1)
    {
      ignore_unused(a0);
      return  boost::simd::splat<result_type>(N);
    }
  };
} } }

#endif
