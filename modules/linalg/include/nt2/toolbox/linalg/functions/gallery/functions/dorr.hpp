/*******************************************************************************
 *         Copyright 2003-2012 LASME UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_FUNCTIONS_DORR_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_FUNCTIONS_DORR_HPP_INCLUDED
#include <nt2/include/functions/eye.hpp>
#include <nt2/include/functions/zeros.hpp>
#include <nt2/include/functions/colon.hpp>
#include <nt2/include/functions/toeplitz.hpp>
#include <nt2/include/functions/pow.hpp>

namespace nt2
{
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::dorr_, tag::cpu_,
                                     (A0)(T), 
                                     (scalar_<integer_<A0> > )
                                     (scalar_<floating_<A1> > )
                                     )
  { 
    typedef T value_type; 
    typedef table< value_type, matlab_index_> result_type; 
    inline result_type operator()(const A0 & n, const A1& theta) const
      {
        result_type c = zeros(n, 1, meta::as_<value_type>);
        result_type d = c;
        result_type e = c; 
        // All length n for convenience.  Make c, e of length n-1 later.
        value_type h = rec(oneplus(value_type(n)));
        size_t m = (n+1)/2; //size_t(nt2::floor ( (n+1)/2.0 ));
        value_type term = theta/nt2::sqr(h);
        
        result_type i = colon(One<value_type>(), value_type(m));
        c(i) = -term;
        e(i) = c(i) - (Half<value_type>()-i*h)/h;
        d(i) = -(c(i) + e(i));
        
        i = colon(value_type(m+1), value_type(n));
        e(i) = -term;
        c(i) = e(i) + (Half<value_type>()-i*h)/h;
        d(i) = -(c(i) + e(i));
        
        return  tridiag (c(AllFrom(2)), d, e(AllTo(n-1)));
      }
  };
}
    
#endif
