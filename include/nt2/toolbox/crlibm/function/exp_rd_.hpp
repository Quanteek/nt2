//////////////////////////////////////////////////////////////////////////////
///   Copyright 2003 and onward LASMEA UMR 6602 CNRS/U.B.P Clermont-Ferrand
///   Copyright 2009 and onward LRI    UMR 8623 CNRS/Univ Paris Sud XI
///
///          Distributed under the Boost Software License, Version 1.0
///                 See accompanying file LICENSE.txt or copy at
///                     http://www.boost.org/LICENSE_1_0.txt
//////////////////////////////////////////////////////////////////////////////
#ifndef NT2_TOOLBOX_CRLIBM_FUNCTION_EXP_RD_HPP_INCLUDED
#define NT2_TOOLBOX_CRLIBM_FUNCTION_EXP_RD_HPP_INCLUDED
#include <nt2/include/simd.hpp>
#include <nt2/include/functor.hpp>
#include <nt2/toolbox/crlibm/include.hpp>

namespace nt2 { namespace crlibm
  {	       
    struct exp_rd__ {};
    NT2_FUNCTION_IMPLEMENTATION(crlibm::exp_rd__, exp_rd_, 1)
  }
}
 
#include <nt2/toolbox/crlibm/function/scalar/exp_rd_.hpp>
// #include NT2_CRLIBM_INCLUDE(exp_rd_.hpp) 

 
#endif
