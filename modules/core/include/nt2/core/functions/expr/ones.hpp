//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#ifndef NT2_CORE_FUNCTIONS_EXPR_ONES_HPP_INCLUDED
#define NT2_CORE_FUNCTIONS_EXPR_ONES_HPP_INCLUDED

#include <nt2/core/container/dsl.hpp>
#include <nt2/core/functions/ones.hpp>
#include <nt2/include/functions/box.hpp>
#include <nt2/core/functions/of_size.hpp>
#include <nt2/include/functions/isrow.hpp>
#include <nt2/include/functions/length.hpp>

namespace nt2 { namespace ext
{
  //============================================================================
  // Generates ones from expression (support size(a) + type calls)
  //============================================================================
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::ones_, tag::cpu_
                            , (A0)
                            , (ast_<A0>)
                            )
  {
    typedef typename  boost::proto::
                      result_of::make_expr< nt2::tag::ones_
                                          , container::domain
                                          , box<of_size_max>
                                          , box< meta::constant_<nt2::tag::One> >
                                          , meta::as_<double>
                                          >::type             result_type;

    BOOST_FORCEINLINE result_type operator()(A0 const& a0) const
    {
      // Expression must be a row vector
      BOOST_ASSERT_MSG
      ( nt2::isrow(a0)
      , "Error using ones: Size vector must be a row vector."
      );

      of_size_max sizee;
      std::size_t sz = std::min(of_size_max::size(),nt2::length(a0));
      for(std::size_t i=1;i<sz;++i) sizee[i-1] = a0(i);

      return boost::proto::make_expr< nt2::tag::ones_
                                    , container::domain
                                    > ( boxify(sizee)
                                      , boxify(meta::constant_<nt2::tag::One>())
                                      , meta::as_<double>()
                                      );
    }
  };

  //============================================================================
  // Generates ones from expression (support size(a) + type calls)
  //============================================================================
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::ones_, tag::cpu_
                            , (A0)(T)
                            , (ast_<A0>)
                              (target_< scalar_< unspecified_<T> > >)
                            )
  {
    typedef typename  boost::proto::
                      result_of::make_expr< nt2::tag::ones_
                                          , container::domain
                                          , box<of_size_max>
                                          , box< meta::constant_<nt2::tag::One> >
                                          , T
                                          >::type             result_type;

    BOOST_FORCEINLINE result_type operator()(A0 const& a0, T const& tgt) const
    {
      // Expression must be a row vector
      BOOST_ASSERT_MSG
      ( nt2::isrow(a0)
      , "Error using ones: Size vector must be a row vector."
      );

      of_size_max sizee;
      std::size_t sz = std::min(of_size_max::size(),nt2::length(a0));
      for(std::size_t i=1;i<sz;++i) sizee[i-1] = a0(i);

      return boost::proto::make_expr< nt2::tag::ones_
                                    , container::domain
                                    > ( boxify(sizee)
                                      , boxify(meta::constant_<nt2::tag::One>())
                                      , tgt
                                      );
    }
  };
} }

#endif
