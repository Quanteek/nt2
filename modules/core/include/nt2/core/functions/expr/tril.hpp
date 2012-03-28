//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#ifndef NT2_CORE_FUNCTIONS_EXPR_TRIL_HPP_INCLUDED
#define NT2_CORE_FUNCTIONS_EXPR_TRIL_HPP_INCLUDED

#include <nt2/core/container/dsl.hpp>
#include <nt2/core/functions/tril.hpp>
#include <nt2/include/functions/box.hpp>
#include <nt2/include/functions/ismatrix.hpp>
#include <nt2/core/functions/details/tril.hpp>

namespace nt2 { namespace ext
{
  //============================================================================
  // Unary tril
  //============================================================================
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::tril_, tag::cpu_, (A0), (ast_<A0>) )
  {
    typedef typename  boost::proto::
                      result_of::make_expr< nt2::tag::tril_
                                          , container::domain
                                          , A0 const&
                                          , box<nt2::details::tril>
                                          >::type             result_type;

    BOOST_FORCEINLINE result_type operator()(A0 const& a0) const
    {
      // Expression must be a matrix
      BOOST_ASSERT_MSG( nt2::ismatrix(a0)
                      , "Error using tril: First input must be 2D."
                      );

      return boost::proto::make_expr< nt2::tag::tril_
                                    , container::domain
                                    > ( boost::cref(a0)
                                      , boxify(nt2::details::tril())
                                      );
    }
  };

  //============================================================================
  // Binary tril
  //============================================================================
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::offset_tril_, tag::cpu_
                            , (A0)(A1)
                            , (ast_<A0>)
                              (scalar_< integer_<A1> >)
                            )
  {
    typedef typename  boost::proto::
                      result_of::make_expr< nt2::tag::offset_tril_
                                          , container::domain
                                          , A0 const&
                                          , box<nt2::details::offset_tril>
                                          >::type             result_type;

    BOOST_FORCEINLINE result_type operator()(A0 const& a0, A1 const& a1) const
    {
      // Expression must be a matrix
      BOOST_ASSERT_MSG( nt2::ismatrix(a0)
                      , "Error using tril: First input must be 2D."
                      );

      return boost::proto::make_expr< nt2::tag::offset_tril_
                                    , container::domain
                                    > ( boost::cref(a0)
                                      , boxify(nt2::details::offset_tril(a1))
                                      );
    }
  };
} }

#endif
