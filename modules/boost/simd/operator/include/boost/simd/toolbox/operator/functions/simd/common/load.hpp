//==============================================================================
//         Copyright 2003 - 2011 LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2011 LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#ifndef BOOST_SIMD_TOOLBOX_OPERATOR_FUNCTIONS_SIMD_COMMON_LOAD_HPP_INCLUDED
#define BOOST_SIMD_TOOLBOX_OPERATOR_FUNCTIONS_SIMD_COMMON_LOAD_HPP_INCLUDED

#include <boost/simd/toolbox/operator/functions/load.hpp>
#include <boost/simd/include/functions/unaligned_load.hpp>
#include <boost/simd/sdk/simd/logical.hpp>
#include <boost/simd/sdk/memory/is_aligned.hpp>
#include <boost/mpl/equal_to.hpp>
#include <boost/assert.hpp>

namespace boost { namespace simd { namespace ext
{
  // regular load
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( boost::simd::tag::load_, tag::cpu_
                            , (A0)(A1)(A2)(X)
                            , (iterator_<scalar_< fundamental_<A0> > >)
                              (scalar_< fundamental_<A1> >)
                              ((target_< simd_< fundamental_<A2>, X > >))
                            )
  {
    typedef typename dispatch::meta::
            call<tag::unaligned_load_(A0 const&, A1 const&, A2 const&)>::type
    result_type;
    inline result_type operator()(const A0& a0, const A1& a1, const A2&) const
    {
      BOOST_ASSERT_MSG
      ( boost::simd::memory::is_aligned(a0,BOOST_SIMD_CONFIG_ALIGNMENT)
     && boost::simd::memory::is_aligned(a0+a1,BOOST_SIMD_CONFIG_ALIGNMENT)
      , "Unaligned memory location. You tried to load with a pointer that"
        " is not aligned on the simd vector size.");
      return unaligned_load<typename A2::type>(a0, a1);
    }
  };

  // shifted load
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( boost::simd::tag::load_, tag::cpu_
                            , (A0)(A1)(A2)(A3)(X)
                            , (iterator_<scalar_< fundamental_<A0> > >)
                              (scalar_< fundamental_<A1> >)
                              ((target_< simd_< fundamental_<A2>, X > >))
                              (mpl_integral_< scalar_< integer_<A3> > >)
                            )
  {
    typedef typename dispatch::meta::
            call<tag::unaligned_load_(A0 const&, A1 const&, A2 const&, A3 const&)>::type
    result_type;
    inline result_type operator()(const A0& a0, const A1& a1,
                                  const A2&, const A3&) const
    {
      BOOST_ASSERT_MSG
      ( boost::simd::memory::is_aligned(a0,BOOST_SIMD_CONFIG_ALIGNMENT)
     && boost::simd::memory::is_aligned(a0+a1,BOOST_SIMD_CONFIG_ALIGNMENT)
      , "Unaligned memory location. You tried to load with a pointer that"
        "is not aligned on the simd vector size.");
      return unaligned_load<typename A2::type, A3::value>(a0, a1);
    }
  };

  // gather
  //TODO Why not a proper gather functor ?
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION_IF( boost::simd::tag::load_, tag::cpu_
                            , (A0)(A1)(A2)(X)
                            , (mpl::equal_to< boost::simd::meta::cardinal_of<A1>
                                            , boost::simd::meta::cardinal_of<typename A2::type>
                                            >
                              )
                            , (iterator_< scalar_< fundamental_<A0> > >)
                              ((simd_< integer_<A1>, X >))
                              ((target_< simd_< fundamental_<A2>, X > >))
                            )
  {
    typedef typename dispatch::meta::
            call<tag::unaligned_load_(A0 const&, A1 const&, A2 const&)>::type
    result_type;
    inline result_type operator()(const A0& a0, const A1& a1, const A2&) const
    {
      BOOST_ASSERT_MSG
      ( boost::simd::memory::is_aligned(a0,BOOST_SIMD_CONFIG_ALIGNMENT)
      , "Unaligned memory location. You tried to load with a pointer that"
        "is not aligned on the simd vector size.");
      return unaligned_load<typename A2::type>(a0, a1);
    }
  };

  // logical
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( boost::simd::tag::load_, tag::cpu_
                            , (A0)(A1)(A2)(X)
                            , (iterator_<scalar_< logical_<A0> > >)
                              (scalar_< fundamental_<A1> >)
                              ((target_< simd_< logical_<A2>, X > >))
                            )
  {
    typedef typename A2::type result_type;
    inline result_type operator()(const A0& a0, const A1& a1, const A2&)const
    {
      BOOST_ASSERT_MSG
      ( boost::simd::memory::is_aligned(a0,BOOST_SIMD_CONFIG_ALIGNMENT)
      , "Unaligned memory location. You tried to load with a pointer that"
        "is not aligned on the simd vector size.");
      return unaligned_load<typename A2::type>(a0, a1);
    }
  };

} } }

#endif
