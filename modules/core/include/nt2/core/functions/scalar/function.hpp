//==============================================================================
//         Copyright 2003 - 2012   LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2012   LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#ifndef NT2_CORE_FUNCTIONS_SCALAR_FUNCTION_HPP_INCLUDED
#define NT2_CORE_FUNCTIONS_SCALAR_FUNCTION_HPP_INCLUDED

#include <nt2/core/functions/function.hpp>
#include <nt2/sdk/parameters.hpp>
#include <nt2/dsl/functions/terminal.hpp>
#include <boost/type_traits/is_const.hpp>
#include <nt2/include/functions/sub2ind.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <boost/fusion/include/vector_tie.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>
#include <boost/preprocessor/repetition/enum_binary_params.hpp>

namespace nt2 { namespace ext
{
#define M0(z,n,t) (I##n)
#define M1(z,n,t) (scalar_< integer_<I##n> >)

#define M2(z,n,t)                                                               \
NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::function_, tag::cpu_                      \
                          , (A0)BOOST_PP_REPEAT(n,M0,~)                         \
                          , (ast_<A0>)BOOST_PP_REPEAT(n,M1,~)                   \
                          )                                                     \
{                                                                               \
  typedef typename boost::remove_const<A0>::type                  base;         \
  typedef typename meta::                                                       \
          scalar_of < typename boost::dispatch::meta::semantic_of<A0&>::type    \
                    >::type                                       result_type;  \
  typedef typename base::index_type::type                         idx_t;        \
                                                                                \
  BOOST_FORCEINLINE result_type                                                 \
  operator()(A0& a0, BOOST_PP_ENUM_BINARY_PARAMS(n,I,i) ) const                 \
  {                                                                             \
    return nt2::run ( a0                                                        \
                    , nt2::sub2ind( nt2::extent(a0)                             \
                          , boost::fusion::                                     \
                            vector_tie(BOOST_PP_ENUM_PARAMS(n,i))               \
                          , idx_t()                                             \
                          )                                                     \
                    , meta::as_<typename base::value_type>()                    \
                    );                                                          \
  }                                                                             \
};                                                                              \
/**/

BOOST_PP_REPEAT_FROM_TO(1,BOOST_PP_INC(NT2_MAX_DIMENSIONS),M2,~)

#undef M0
#undef M1
#undef M2
} }

#endif
