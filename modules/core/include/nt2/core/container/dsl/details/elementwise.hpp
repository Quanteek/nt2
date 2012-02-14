//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#ifndef NT2_CORE_CONTAINER_DSL_DETAILS_ELEMENTWISE_HPP_INCLUDED
#define NT2_CORE_CONTAINER_DSL_DETAILS_ELEMENTWISE_HPP_INCLUDED

#include <nt2/core/container/dsl/generator.hpp>
#include <boost/assert.hpp>
#include <boost/fusion/include/transform.hpp>
#include <boost/fusion/include/fold.hpp>
#include <boost/fusion/include/at_c.hpp>

namespace nt2 { namespace container { namespace ext
{
  // element-wise size selection logic
  struct size_fold
  {
    template<class A0, class A1, class Dummy = void>
    struct select;

    template<class Sig>
    struct result;

    template<class This, class A0, class A1>
    struct result<This(A0, A1)>
    {
      typedef typename
      select< typename meta::strip<A0>::type
            , typename meta::strip<A1>::type
            >::type const& type;
    };

    template<class A1, class Dummy>
    struct select<_0D, A1, Dummy>
    {
      typedef A1 type;
    };
    template<class A1>
    BOOST_DISPATCH_FORCE_INLINE
    typename result<size_fold(_0D const&, A1 const&)>::type
    operator()(_0D const&, A1 const& a1) const
    {
      return a1;
    }

    template<class A0, class Dummy>
    struct select<A0, _0D, Dummy>
    {
      typedef A0 type;
    };
    template<class A0>
    BOOST_DISPATCH_FORCE_INLINE
    typename result<size_fold(A0 const&, _0D const&)>::type
    operator()(A0 const& a0, _0D const&) const
    {
      return a0;
    }

    template<class Dummy>
    struct select<_0D, _0D, Dummy>
    {
      typedef _0D type;
    };
    BOOST_DISPATCH_FORCE_INLINE
    result<size_fold(_0D const&, _0D const&)>::type
    operator()(_0D const& a0, _0D const&) const
    {
      return a0;
    }

    template<class A0, class A1, class Dummy>
    struct select
    {
      typedef typename
      boost::mpl::if_c< A0::static_status,A0,A1>::type type;
    };

    template<class A0, class A1> BOOST_DISPATCH_FORCE_INLINE
    typename result<size_fold(A0 const&, A1 const&)>::type
    operator()(A0 const& a0, A1 const& a1) const
    {
      BOOST_ASSERT_MSG(a0 == a1, "Sizes are not compatible");
      return selection(a0,a1,boost::mpl::bool_<A0::static_status>());
    }

    template<class A0, class A1> BOOST_DISPATCH_FORCE_INLINE
    typename result<size_fold(A0 const&, A1 const&)>::type
    selection(A0 const& a0, A1 const&, boost::mpl::true_ const&) const
    {
      return a0;
    }

    template<class A0, class A1> BOOST_DISPATCH_FORCE_INLINE
    typename result<size_fold(A0 const&, A1 const&)>::type
    selection(A0 const&, A1 const& a1, boost::mpl::false_ const&) const
    {
      return a1;
    }
  };

  // element-wise size n-ary
  template<class Tag, class Domain, int N, class Expr>
  struct size_of
  {
    typedef typename boost::fusion::result_of::
    transform<Expr const, size_transform<Domain> >::type sizes;

    typedef typename boost::fusion::result_of::
    at_c<sizes, 0>::type init;

    typedef typename boost::fusion::result_of::
    fold<sizes, init const, size_fold>::type  result_type;

    BOOST_FORCEINLINE
    result_type operator()(Expr& e) const
    {
      sizes sz = boost::fusion::transform(e, size_transform<Domain>());
      return boost::fusion::fold(sz, boost::fusion::at_c<0>(sz), size_fold());
    }
  };

  // element-wise size unary
  template<class Tag, class Domain, class Expr>
  struct size_of<Tag, Domain, 1, Expr>
  {
    typedef typename boost::proto::result_of::
    child_c<Expr&, 0>::type                         child0;

    typedef typename size_transform<Domain>::template
    result<size_transform<Domain>(child0)>::type    result_type;

    BOOST_FORCEINLINE
    result_type operator()(Expr& e) const
    {
      return size_transform<Domain>()(boost::proto::child_c<0>(e));
    }
  };

} } }

#endif
