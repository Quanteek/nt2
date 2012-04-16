//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#ifndef NT2_DSL_FUNCTIONS_CONTAINER_TERMINAL_HPP_INCLUDED
#define NT2_DSL_FUNCTIONS_CONTAINER_TERMINAL_HPP_INCLUDED

#include <nt2/sdk/simd/category.hpp>
#include <boost/fusion/include/size.hpp>
#include <nt2/dsl/functions/terminal.hpp>
#include <nt2/include/functions/splat.hpp>
#include <nt2/include/functions/unaligned_load.hpp>
#include <nt2/include/functions/unaligned_store.hpp>
#include <nt2/core/settings/details/fusion.hpp>
#include <nt2/core/container/category.hpp>

namespace nt2 { namespace ext
{
  //============================================================================
  // table terminal with a position in scalar read mode
  //============================================================================
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::terminal_, tag::cpu_
                            , (A0)(T0)(S0)(State)(Data)
                            , ((expr_< table_< unspecified_<A0>, S0 >
                                     , T0
                                     , boost::mpl::long_<0>
                                     >
                              ))
                              (generic_< integer_<State> >)
                              (target_<scalar_<unspecified_<Data> > >)
                            )
  {
    typedef typename boost::dispatch::meta::
    scalar_of<A0&>::type                               result_type;

    BOOST_FORCEINLINE result_type
    operator()(A0& a0, State const& state, Data const&) const
    {
       return nt2::terminal(a0)[state];
    }
  };

  //============================================================================
  // table terminal with a position in scalar write mode
  //============================================================================
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::terminal_, tag::cpu_
                            , (A0)(T0)(S0)(State)(Data)
                            , ((expr_< table_< unspecified_<A0>, S0 >
                                     , T0
                                     , boost::mpl::long_<0>
                                     >
                              ))
                              (generic_< integer_<State> >)
                              (scalar_<unspecified_<Data> >)
                            )
  {
    typedef typename boost::dispatch::meta::
    scalar_of<A0&>::type                                result_type;

    BOOST_FORCEINLINE result_type
    operator()(A0& a0, State const& state, Data const& data) const
    {
       return nt2::terminal(a0)[state] = data;
    }
  };

  //============================================================================
  // table terminal with a position in SIMD read mode
  //============================================================================
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::terminal_, tag::cpu_
                            , (A0)(T0)(S0)(State)(Data)(X)
                            , ((expr_< table_< unspecified_<A0>, S0 >
                                     , T0
                                     , boost::mpl::long_<0>
                                     >
                              ))
                              (generic_< integer_<State> >)
                              ((target_< simd_<unspecified_<Data>, X> >))
                            )
  {
    typedef typename boost::dispatch::meta::
            strip< typename boost::dispatch::meta::
                   scalar_of< typename boost::dispatch::meta::
                              semantic_of<A0&>::type
                            >::type
                 >::type                            stype;

    typedef boost::simd::native<stype, X>           result_type;

    BOOST_FORCEINLINE
    result_type operator()(A0 const& a0, State const& state, Data const&) const
    {
      return unaligned_load<result_type>(nt2::terminal(a0).raw(), state);
    }
  };

  //============================================================================
  // table terminal with a position in SIMD write mode
  //============================================================================
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::terminal_, tag::cpu_
                            , (A0)(T0)(S0)(State)(Data)(X)
                            , ((expr_< table_< unspecified_<A0>, S0 >
                                     , T0
                                     , boost::mpl::long_<0>
                                     >
                              ))
                              (generic_< integer_<State> >)
                              ((simd_<unspecified_<Data>, X>))
                            )
  {
    typedef Data                                            result_type;

    BOOST_FORCEINLINE
    result_type operator()(A0& a0, State const& state, Data const& data) const
    {
      return unaligned_store<result_type>(data, nt2::terminal(a0).raw(), state);
    }
  };

  //============================================================================
  // scalar terminal, return value in scalar mode (LHS not allowed)
  //============================================================================
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::terminal_, tag::cpu_
                            , (A0)(T0)(State)(Data)
                            , ((expr_< scalar_< unspecified_<A0> >
                                     , T0
                                     , boost::mpl::long_<0>
                                     >
                              ))
                              (generic_< integer_<State> >)
                              (target_< scalar_< unspecified_<Data> > >)
                            )
  {
    typedef typename nt2::meta::call<T0(A0&)>::type result_type;

    BOOST_FORCEINLINE result_type
    operator()(A0& a0, State const&, Data const&) const
    {
      return nt2::terminal(a0);
    }
  };

  //============================================================================
  // scalar terminal, splat value in SIMD read mode (LHS not allowed)
  //============================================================================
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::terminal_, tag::cpu_
                            , (A0)(T0)(State)(Data)(X)
                            , ((expr_< scalar_< unspecified_<A0> >
                                     , T0
                                     , boost::mpl::long_<0>
                                     >
                              ))
                              (generic_< integer_<State> >)
                              ((target_< simd_< unspecified_<Data>,X > >))
                            )
  {
    typedef typename boost::dispatch::meta::
            strip< typename boost::dispatch::meta::
                   semantic_of<A0&>::type
                 >::type                            stype;

    typedef boost::simd::native<stype, X>           result_type;

    BOOST_FORCEINLINE
    result_type operator()(A0& a0, State const&, Data const&) const
    {
      return nt2::splat<result_type>(nt2::terminal(a0));
    }
  };
} }

#endif
