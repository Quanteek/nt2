//==============================================================================
//         Copyright 2003 - 2011 LASMEA UMR 6602 CNRS/Univ. Clermont II         
//         Copyright 2009 - 2011 LRI    UMR 8623 CNRS/Univ Paris Sud XI         
//                                                                              
//          Distributed under the Boost Software License, Version 1.0.          
//                 See accompanying file LICENSE.txt or copy at                 
//                     http://www.boost.org/LICENSE_1_0.txt                     
//==============================================================================
#ifndef BOOST_SIMD_TOOLBOX_OPERATOR_FUNCTIONS_SCALAR_LOAD_HPP_INCLUDED
#define BOOST_SIMD_TOOLBOX_OPERATOR_FUNCTIONS_SCALAR_LOAD_HPP_INCLUDED
#include <boost/dispatch/meta/mpl.hpp>
#include <boost/simd/sdk/memory/details/category.hpp>
#include <boost/dispatch/functor/preprocessor/call.hpp>

namespace boost { namespace simd { namespace ext
{
//==============================================================================
// load_ without offset
//==============================================================================
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( boost::simd::tag::load_ , tag::cpu_
                            , (A0)(A1)(A2)
                            , (iterator_< scalar_< fundamental_<A0> > >)
                              (scalar_< fundamental_<A1> >)
                              (target_< scalar_< fundamental_<A2> > >)
                            )
  {
    typedef typename A2::type result_type;
    inline result_type operator()(const A0& a0, const A1& a1, const A2&) const
    {
      A0 that = a0;
      std::advance(that,a1);
      return *that;
    }
  };

//==============================================================================
// load_ with offset
//==============================================================================
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION(  boost::simd::tag::load_ , tag::cpu_
                            , (A0)(A1)(A2)(A3)
                            , (iterator_< scalar_< fundamental_<A0> > >)
                              (scalar_< fundamental_<A1> >)
                              (target_< scalar_< fundamental_<A2> > >)
                              (mpl_integral_< scalar_< integer_<A3> > >)
                            )
  {
    typedef typename A2::type result_type;
    inline result_type operator()(const A0& a0, const A1& a1,
                                  const A2&, const A3&) const

    {
      A0 that = a0;
      std::advance(that,a1+A3::value);
      return *that;
    }
  };
} } }

#endif
