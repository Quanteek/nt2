//////////////////////////////////////////////////////////////////////////////
///   Copyright 2003 and onward LASMEA UMR 6602 CNRS/U.B.P Clermont-Ferrand
///   Copyright 2009 and onward LRI    UMR 8623 CNRS/Univ Paris Sud XI
///
///          Distributed under the Boost Software License, Version 1.0
///                 See accompanying file LICENSE.txt or copy at
///                     http://www.boost.org/LICENSE_1_0.txt
//////////////////////////////////////////////////////////////////////////////
#ifndef NT2_TOOLBOX_FUZZY_FUNCTIONS_FUZZY_DEFINITELY_LESS_HPP_INCLUDED
#define NT2_TOOLBOX_FUZZY_FUNCTIONS_FUZZY_DEFINITELY_LESS_HPP_INCLUDED
#include <nt2/include/simd.hpp>
#include <nt2/include/functor.hpp>
#include <nt2/toolbox/fuzzy/include.hpp>

namespace nt2 { namespace tag
  {         
    struct fuzzy_definitely_less_ {};
  }
  NT2_FUNCTION_IMPLEMENTATION(tag::fuzzy_definitely_less_, fuzzy_definitely_less, 3)
}
 
 
#endif

// modified by jt the 25/12/2010