//==============================================================================
//         Copyright 2003 - 2011 LASMEA UMR 6602 CNRS/Univ. Clermont II         
//         Copyright 2009 - 2011 LRI    UMR 8623 CNRS/Univ Paris Sud XI         
//                                                                              
//          Distributed under the Boost Software License, Version 1.0.          
//                 See accompanying file LICENSE.txt or copy at                 
//                     http://www.boost.org/LICENSE_1_0.txt                     
//==============================================================================
#ifndef NT2_TOOLBOX_HYPERBOLIC_FUNCTIONS_SCALAR_SINHCOSH_HPP_INCLUDED
#define NT2_TOOLBOX_HYPERBOLIC_FUNCTIONS_SCALAR_SINHCOSH_HPP_INCLUDED
#include <nt2/sdk/meta/as_logical.hpp>
#include <boost/fusion/tuple.hpp>
#include <nt2/include/functions/tofloat.hpp>
#include <nt2/include/functions/expm1.hpp>
#include <nt2/include/functions/if_else.hpp>
#include <nt2/include/functions/is_eqz.hpp>
#include <nt2/include/functions/is_equal.hpp>
#include <nt2/include/functions/oneplus.hpp>
#include <nt2/include/functions/abs.hpp>
#include <nt2/include/constants/inf.hpp>
#include <nt2/include/constants/half.hpp>

namespace nt2 { namespace ext
{
  NT2_FUNCTOR_IMPLEMENTATION(  nt2::tag::sinhcosh_, tag::cpu_,(A0)(A1)
                               , ((scalar_<arithmetic_<A0> >))
                                 ((scalar_<floating_<A1> >))
                                 ((scalar_<floating_<A1> >))
                             )
  {
    typedef int result_type;    
    inline result_type operator()(A0 const& a0,A1 & a1,A1 & a2) const
    {
      if (is_inf(a0))
        {
          a1 = a0;
          a2 = nt2::abs(a0); 
        }
      else
        {
          const A1 u = expm1(tofloat(a0));
          const A1 up1 = oneplus(u);
          const A1 tmp =u/up1; 
          a1 = Half<A1>()*tmp*(oneplus(up1));
          a2 = oneplus(Half<A1>()*tmp*u);
        }
      return 0;
    }
  };
 
  NT2_FUNCTOR_IMPLEMENTATION(nt2::tag::sinhcosh_, tag::cpu_,(A0)(A1),
                             ((scalar_ < arithmetic_<A0> > ))
                             ((scalar_ < floating_<A1> > ))
                             )
    {
    typedef A1 result_type;    
    inline result_type operator()(A0 const& a0,A1 & a2) const
    {
      A1 a1;
      sinhcosh(tofloat(a0),a1, a2);
      return a1;
    }
  };

  NT2_FUNCTOR_IMPLEMENTATION(nt2::tag::sinhcosh_, tag::cpu_,
                         (A0),
                         ((scalar_<arithmetic_<A0> >))
                        )
  {
      typedef typename meta::as_floating<A0>::type  rtype;
      typedef boost::fusion::tuple<rtype, rtype> result_type;
    
    NT2_FUNCTOR_CALL(1)
    {
      result_type res;
      sinhcosh(tofloat(a0), boost::fusion::at_c<0>(res), 
                 boost::fusion::at_c<1>(res));
      return res;
    }
  };
} }
#endif
