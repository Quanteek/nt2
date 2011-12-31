//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#ifndef NT2_TOOLBOX_COMPLEX_CONSTANTS_GENERIC_I_HPP_INCLUDED
#define NT2_TOOLBOX_COMPLEX_CONSTANTS_GENERIC_I_HPP_INCLUDED

#include <nt2/toolbox/complex/constants/i.hpp>
#include <nt2/include/constants/one.hpp>
#include <nt2/sdk/complex/imaginary.hpp>
#include <nt2/sdk/complex/meta/as_imaginary.hpp>

namespace nt2 { namespace ext
{
  NT2_FUNCTOR_IMPLEMENTATION ( nt2::tag::I, tag::cpu_, (A0)
                             , ((target_< generic_< imaginary_< arithmetic_<A0> > > >))
                             )
  {
    typedef typename A0::type base_type;
    //    typedef typename meta::as_imaginary<base_type>::type result_type;
    typedef                   imaginary<base_type>       result_type;      
    BOOST_DISPATCH_FORCE_INLINE result_type operator()(A0 const&) const
    {
      return result_type(One<typename result_type::type>());
    }
  };
} }

#endif