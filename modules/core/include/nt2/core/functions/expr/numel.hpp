//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#ifndef NT2_CORE_FUNCTIONS_EXPR_NUMEL_HPP_INCLUDED
#define NT2_CORE_FUNCTIONS_EXPR_NUMEL_HPP_INCLUDED

#include <nt2/core/functions/numel.hpp>
#include <nt2/include/functions/extent.hpp>
#include <nt2/include/functions/multiplies.hpp>
#include <boost/fusion/include/fold.hpp>
#include <nt2/core/container/dsl.hpp>
#include <nt2/core/container/category.hpp>

namespace nt2 { namespace ext
{
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::numel_, tag::cpu_
                            , (A0)
                            , (fusion_sequence_<A0>)
                            )
  {
    typedef typename boost::fusion::result_of::
    fold< A0
        , boost::mpl::size_t<1>
        , nt2::functor<tag::multiplies_>
        >::type                                                   result_type;

    BOOST_DISPATCH_FORCE_INLINE
    result_type operator()(const A0& a0) const
    {
      return boost::fusion::fold( a0
                                , boost::mpl::size_t<1>()
                                , functor<tag::multiplies_>()
                                );
    }
  };

  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::numel_, tag::cpu_
                            , (A0)(T)
                            , ((expr_ < unspecified_<A0>
                                      , nt2::container::domain, T
                                      >
                              ))
                            )
  {
    typedef typename meta::call<tag::numel_(typename A0::extent_type)>::type result_type;

    BOOST_DISPATCH_FORCE_INLINE
    result_type operator()(const A0& a0) const
    {
      return nt2::numel(nt2::extent(a0));
    }
  };
} }


#endif
