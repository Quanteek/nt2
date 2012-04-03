/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_CORE_FUNCTIONS_SCALED_ROWS_HPP_INCLUDED
#define NT2_CORE_FUNCTIONS_SCALED_ROWS_HPP_INCLUDED
#include <nt2/include/functor.hpp>
#include <nt2/core/functions/details/scaled_rows.hpp>
#include <nt2/sdk/meta/generative_hierarchy.hpp>
#include <nt2/core/container/dsl/details/generative.hpp>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>

namespace nt2
{
  namespace tag
  {
    struct scaled_rows_ : ext::generative_<scaled_rows_>
    {
      typedef ext::generative_<scaled_rows_> parent;
    };
  }

  NT2_FUNCTION_IMPLEMENTATION(nt2::tag::scaled_rows_, scaled_rows, 3)
  NT2_FUNCTION_IMPLEMENTATION(nt2::tag::scaled_rows_, scaled_rows, 4)
}

namespace nt2 { namespace container { namespace ext
{
  //============================================================================
  // Register scaled_rows as a generative expression
  //============================================================================
  template<class Domain, class Expr, int N>
  struct generator<tag::scaled_rows_,Domain,N,Expr>   : generative_generator<Expr>
  {};

  template<class Domain, class Expr, int N>
  struct size_of<tag::scaled_rows_,Domain,N,Expr>     : generative_size_of<Expr>
  {};
} } }

#endif

// /////////////////////////////////////////////////////////////////////////////
// End of scaled_rows.hpp
// /////////////////////////////////////////////////////////////////////////////
