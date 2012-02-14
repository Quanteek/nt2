//==============================================================================
//         Copyright 2003 - 2011 LASMEA UMR 6602 CNRS/Univ. Clermont II         
//         Copyright 2009 - 2011 LRI    UMR 8623 CNRS/Univ Paris Sud XI         
//                                                                              
//          Distributed under the Boost Software License, Version 1.0.          
//                 See accompanying file LICENSE.txt or copy at                 
//                     http://www.boost.org/LICENSE_1_0.txt                     
//==============================================================================
#ifndef NT2_TOOLBOX_HYPERBOLIC_FUNCTIONS_COMPLEX_GENERIC_SINH_HPP_INCLUDED
#define NT2_TOOLBOX_HYPERBOLIC_FUNCTIONS_COMPLEX_GENERIC_SINH_HPP_INCLUDED
#include <nt2/include/functions/sincos.hpp>
#include <nt2/include/functions/sinhcosh.hpp>
#include <nt2/include/functions/sinh.hpp>
#include <nt2/include/functions/cosh.hpp>
#include <nt2/include/functions/real.hpp>
#include <nt2/include/functions/imag.hpp>
#include <nt2/include/functions/sin.hpp>
#include <nt2/include/functions/is_eqz.hpp>
#include <nt2/include/functions/sign.hpp>
#include <nt2/include/functions/abs.hpp>
#include <nt2/sdk/complex/meta/as_complex.hpp>
#include <nt2/sdk/complex/meta/as_real.hpp>
#include <nt2/sdk/complex/meta/as_dry.hpp>
#include <nt2/include/functions/bitwise_cast.hpp>
// #include <iostream>

//sinh(x+iy)=sinh(x)cos(y)+i.cosh(x)sin(y).
namespace nt2 { namespace ext
{
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::sinh_, tag::cpu_, (A0)
                            , (generic_< complex_< arithmetic_<A0> > >)
                            )
  {
    typedef A0 result_type;
    NT2_FUNCTOR_CALL(1)
    {
      typedef typename meta::as_real<A0>::type rtype; 
      rtype c, s, ch, sh;
      sincos(imag(a0), s, c);
      sinhcosh(real(a0), sh, ch);
      rtype r = c*sh;
      rtype i = s*ch;
//       std::cout << "a0   " << a0<< std::endl;                
//       std::cout << "c    " << c << std::endl; 
//       std::cout << "s    " << s << std::endl; 
//       std::cout << "ch   " << ch << std::endl; 
//       std::cout << "sh   " << sh << std::endl; 
//       std::cout << "c*sh = r    " << r << std::endl; 
//       std::cout << "s*ch = i    " << i << std::endl; 
      if (none(is_invalid(a0))) return result_type(r, i);
      r = if_else(logical_and(is_inf(real(a0)), is_invalid(imag(a0))), real(a0), r);
      i = if_else(logical_and(is_inf(real(a0)), is_nan(imag(a0))), nt2::Nan<rtype>(), i);
      r = if_else(is_nan(real(a0)), real(a0), r);
      i = if_else(is_nan(real(a0)), real(a0), i);
      i = if_zero_else(is_real(a0), i);
      r = if_zero_else(is_imag(a0), r);
      result_type res =  result_type(r, i);
      return res;

    }
  };

  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::sinh_, tag::cpu_, (A0)
                            , (generic_< imaginary_< arithmetic_<A0> > >)
                            )
  {
    typedef A0 result_type; 
    NT2_FUNCTOR_CALL(1)
    {
      return bitwise_cast<result_type>(nt2::sin(imag(a0))); 
    }
  };

  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::sinh_, tag::cpu_, (A0)
                            , (generic_< dry_< arithmetic_<A0> > >)
                            )
  {
    typedef A0 result_type; 
    NT2_FUNCTOR_CALL(1)
    {
      return bitwise_cast<result_type>(nt2::sinh(real(a0))); 
    }
  };
  
} }

#endif
