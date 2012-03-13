//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#ifndef NT2_CORE_FUNCTIONS_EYE_HPP_INCLUDED
#define NT2_CORE_FUNCTIONS_EYE_HPP_INCLUDED

#include <nt2/include/functor.hpp>
#include <nt2/core/functions/details/eye.hpp>
#include <nt2/sdk/meta/generative_hierarchy.hpp>
#include <nt2/core/container/dsl/details/generative.hpp>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>

namespace nt2
{
  namespace tag
  {
    struct eye_ : ext::generative_<eye_>
    {
      typedef ext::generative_<eye_> parent;
    };
  }

  NT2_FUNCTION_IMPLEMENTATION(nt2::tag::eye_, eye, 1)
  NT2_FUNCTION_IMPLEMENTATION(nt2::tag::eye_, eye, 2)
  NT2_FUNCTION_IMPLEMENTATION(nt2::tag::eye_, eye, 3)
}

namespace nt2 { namespace container { namespace ext
{
  //============================================================================
  // Register eye as a generative expression
  //============================================================================
  template<class Domain, class Expr, int N>
  struct generator<tag::eye_,Domain,N,Expr>   : generative_generator<Expr>
  {};

  template<class Domain, class Expr, int N>
  struct size_of<tag::eye_,Domain,N,Expr>     : generative_size_of<Expr>
  {};
} } }

#endif
