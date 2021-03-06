//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#ifndef NT2_CORE_CONTAINER_DSL_DETAILS_REDUCTION_HPP_INCLUDED
#define NT2_CORE_CONTAINER_DSL_DETAILS_REDUCTION_HPP_INCLUDED

#include <nt2/core/container/dsl/forward.hpp>
#include <nt2/core/utility/of_size/predef.hpp>
#include <nt2/sdk/memory/container.hpp>
#include <nt2/sdk/meta/strip.hpp>
#include <boost/proto/traits.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <nt2/include/functions/firstnonsingleton.hpp>
#include <nt2/include/functions/terminal.hpp>

namespace nt2 { namespace container { namespace ext
{
  //============================================================================
  // This is the factorized size_of for all reduction function.
  // For any given reduction function tag RED, the registration of their
  // size_of is simply :
  //
  // namespace nt2 { namespace container { namespace ext
  // {
  //  template<class Domain, class Expr, int N>
  //  struct size_of<RED,Domain,N,Expr> : reduction_size_of<Expr, N>
  //  {};
  // } } }
  //
  //============================================================================
  template<class RED, int N, class Expr> struct reduction_size_of;

  template<class RED, class Expr> struct reduction_size_of<RED,2,Expr>
  {
    typedef typename boost::proto::result_of::child_c<Expr&, 0>::type    child0;
    typedef typename meta::strip<child0>::type                          schild0;
    typedef typename meta::strip<typename schild0::extent_type>::type     ext_t;
    typedef typename nt2::make_size<ext_t::static_size>::type       result_type;

    BOOST_FORCEINLINE result_type operator()(Expr& e) const
    {
      result_type res = boost::proto::child_c<0>(e).extent();
      std::size_t red_dim = nt2::terminal(boost::proto::child_c<1>(e)) - 1;

      // If we reduce over the number of dimensions, do nothing
      if(red_dim < result_type::static_size) res[red_dim] = 1;

      return res;
    }
  };

  template<class RED, class Expr> struct reduction_size_of<RED,1,Expr>
  {
    typedef typename boost::proto::result_of::child_c<Expr&, 0>::type child0;
    typedef typename meta::strip<child0>::type                        schild0;
    typedef typename meta::strip<typename schild0::extent_type>::type result_type;

    BOOST_FORCEINLINE result_type operator()(Expr& e) const
    {
      result_type res = boost::proto::child_c<0>(e).extent();
      res[nt2::firstnonsingleton(res)-1] = 1;
      return res;
    }
  };

  //============================================================================
  // This is the factorized generator for all reduction function.
  // For any given reduction function tag RED, the registration of their
  // reshaper is simply :
  //
  // namespace nt2 { namespace container { namespace ext
  // {
  //  template<class Domain, class Expr, int N>
  //  struct generator<RED,Domain,N,Expr> : reduction_generator<Expr>
  //  {};
  // } } }
  //
  //============================================================================
  template<class RED, int N, class Expr> struct reduction_generator
  {
    typedef typename boost::proto::result_of::child_c<Expr&, 0>::type             expr_t;
    typedef typename boost::dispatch::meta::semantic_of<expr_t>::type             sema_t;
    typedef typename meta::strip<sema_t>::type                                    ssema_t;

    typedef typename reduction_size_of<RED,N,Expr>::result_type                   size_type;
    typedef typename nt2::memory::container< typename ssema_t::value_type
                                           , size_type
                                           >                                      type;

    typedef expression< typename boost::remove_const<Expr>::type
                      , type
                      >                                                           result_type;

    BOOST_FORCEINLINE result_type operator()(Expr& e) const
    {
      return result_type(e);
    }
  };
} } }

#endif
