//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#ifndef NT2_CORE_FUNCTIONS_COMMON_AREOFSAMEDIMENSIONS_HPP_INCLUDED
#define NT2_CORE_FUNCTIONS_COMMON_AREOFSAMEDIMENSIONS_HPP_INCLUDED

#include <nt2/core/functions/areofsamedimensions.hpp>
#include <nt2/include/functions/ismatrix.hpp>

namespace nt2 { namespace ext
{
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::areofsamedimensions_, tag::cpu_
                              , (A0)(A1)
                              , (unspecified_<A0>)
                              (unspecified_<A1>)
                              )
  {
    typedef bool result_type;

    BOOST_DISPATCH_FORCE_INLINE result_type operator()(const A0& a0, const A1& a1) const
    {
      typename meta::call<tag::extent_(A0 const&)>::type ex0 = nt2::extent(a0);
      typename meta::call<tag::extent_(A1 const&)>::type ex1 = nt2::extent(a1);

      return ex0 == ex1;
    }
  };
} }

#endif
