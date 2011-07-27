//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#ifndef BOOST_SIMD_TOOLBOX_CONSTANT_CONSTANTS_DETAILS_IEEE_SPEC_HPP_INCLUDED
#define BOOST_SIMD_TOOLBOX_CONSTANT_CONSTANTS_DETAILS_IEEE_SPEC_HPP_INCLUDED

#include <boost/dispatch/meta/strip.hpp>
#include <boost/dispatch/meta/as_integer.hpp>
#include <boost/simd/include/functions/splat.hpp>
#include <boost/dispatch/functor/preprocessor/call.hpp>

#define LOCAL_CONST(TAG, D, F)                                        \
BOOST_DISPATCH_FUNCTOR_IMPLEMENTATION( TAG,tag::cpu_,(A0)                        \
                          , (target_< scalar_< double_<A0> > > )      \
                          )                                           \
{                                                                     \
  typedef typename as_integer < typename A0::type                     \
                              , signed                                \
                              >::type result_type;                    \
  BOOST_DISPATCH_FUNCTOR_CALL(1)                                                 \
  {                                                                   \
    ignore_unused(a0);                                                \
    return boost::simd::splat<result_type>(D);                                     \
  }                                                                   \
};                                                                    \
                                                                      \
BOOST_DISPATCH_FUNCTOR_IMPLEMENTATION( TAG,tag::cpu_,(A0)                        \
                          , (target_< scalar_< float_<A0> > > )       \
                          )                                           \
{                                                                     \
  typedef typename as_integer < typename A0::type                     \
                              , signed                                \
                              >::type result_type;                    \
  BOOST_DISPATCH_FUNCTOR_CALL(1)                                                 \
  {                                                                   \
    ignore_unused(a0);                                                \
    return boost::simd::splat<result_type>(F);                                     \
  }                                                                   \
};                                                                    \
/**/

namespace boost { namespace dispatch { namespace meta
{
  LOCAL_CONST(boost::simd::tag::nb_mantissa_bits_,                  52,         23);
  LOCAL_CONST(boost::simd::tag::nb_exponent_bits_,                  11,          8);
  LOCAL_CONST(boost::simd::tag::max_exponent_    ,                1023,        127);
  LOCAL_CONST(boost::simd::tag::min_exponent_    ,               -1022,       -126);
  LOCAL_CONST(boost::simd::tag::nb_digits_       ,                  53,         24);
  LOCAL_CONST(boost::simd::tag::ldexp_mask_      ,0x7FF0000000000000ll, 0x7F800000);
} } }

#undef LOCAL_CONST

#endif
