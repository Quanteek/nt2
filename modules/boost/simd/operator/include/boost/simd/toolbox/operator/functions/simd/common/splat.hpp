//==============================================================================
//         Copyright 2003 - 2011 LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2011 LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#ifndef BOOST_SIMD_TOOLBOX_OPERATOR_FUNCTIONS_SIMD_COMMON_SPLAT_HPP_INCLUDED
#define BOOST_SIMD_TOOLBOX_OPERATOR_FUNCTIONS_SIMD_COMMON_SPLAT_HPP_INCLUDED

#include <boost/simd/toolbox/operator/functions/splat.hpp>
#include <boost/simd/include/functions/insert.hpp>
#include <boost/simd/sdk/meta/cardinal_of.hpp>
#include <boost/simd/sdk/meta/scalar_of.hpp>

namespace boost { namespace simd { namespace ext
{
  //============================================================================
  // With no idea what we're doing, just fill the vector piecewise
  //============================================================================
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( boost::simd::tag::splat_, tag::cpu_
                            , (A0)(A1)(X)
                            , (unspecified_<A0>)
                              ((target_< simd_< arithmetic_<A1>, X > >))
                            )
  {
    typedef typename A1::type result_type;

    inline result_type operator()(const A0& a0, const A1&) const
    {
      typedef typename meta::scalar_of<result_type>::type sA1;

      result_type tmp;
      for(unsigned int i = 0; i != boost::simd::meta::cardinal_of<result_type>::value; ++i)
        tmp[i] = static_cast<sA1>(a0);

      return tmp;
    }
  };

  //============================================================================
  // Splatting a SIMD value is identity
  //============================================================================
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( boost::simd::tag::splat_, tag::cpu_
                            , (A0)(A1)(X)
                            , ((simd_< arithmetic_<A0>, X >))
                              ((target_< simd_< arithmetic_<A1>, X > >))
                            )
  {
    typedef A0 result_type;

    inline result_type operator()(const A0& a0, const A1&) const
    {
      return a0;
    }
  };
} } }

#endif
