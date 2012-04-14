/*******************************************************************************
 *         Copyright 2003-2012 LASME UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_FUNCTIONS_CHOW_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_FUNCTIONS_CHOW_HPP_INCLUDED
#include <nt2/include/functions/eye.hpp>
#include <nt2/include/functions/zeros.hpp>
#include <nt2/include/functions/colon.hpp>
#include <nt2/include/functions/toeplitz.hpp>
#include <nt2/include/functions/pow.hpp>

namespace nt2
{
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::chow_, tag::cpu_,
                                     (A0)(T), 
                                     (scalar_<integer_<A0> > )
                                     (target_<scalar_<floating_<T> > >)
                                     )
  { 
    typedef typename T::value_type value_type; 
    typedef table< value_type > result_type; 
    inline result_type operator()(const A0 & n, const T&) const
      {
        return chow(n, One<value_type>(), Zero<value_type>()); 
      }
  };
  
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::chow_, tag::cpu_,
                                     (A0)(A1)(A2),
                                     (scalar_<integer_ <A0> )
                                     (scalar_<floating_ <A1> )
                                     (scalar_<floating_ <A2> )
                                     )
  {
    typedef table<A1> result_type; 
    inline result_type operator()(const A0 & n, const A1& alpha,
                                  const A2 & delta) const
      { 
        typedef meta::as_integer<A1>::type i_type; 
        result_type z2 = zeros(1, n); z2(0) = alpha;  z2(1) = 1; 
        result_type z1 = pow(alpha, _(itype(1), itype(n))); 
        result_type r =  toeplitz( z1, z2 ) + delta*eye(n, n);
        return r; 
      };  
  }
}

#endif
