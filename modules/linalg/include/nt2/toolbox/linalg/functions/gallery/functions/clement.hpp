/*******************************************************************************
 *         Copyright 2003-2012 LASME UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_FUNCTIONS_CLEMENT_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_FUNCTIONS_CLEMENT_HPP_INCLUDED
#include <nt2/include/functions/eye.hpp>
#include <nt2/include/functions/zeros.hpp>
#include <nt2/include/functions/colon.hpp>
#include <nt2/include/functions/toeplitz.hpp>
#include <nt2/include/functions/pow.hpp>

namespace nt2
{
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::clement_, tag::cpu_,
                                     (A0)(T), 
                                     (scalar_<integer_<A0> > )
                                     (target_<scalar_<floating_<T> > >)
                                     )
  { 
    typedef typename T::value_type value_type; 
    typedef table< value_type > result_type; 
    inline result_type operator()(const A0 & n, const T&) const
      {
        return clement(n, 0, T()); 
      }
  };
  
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::clement_, tag::cpu_,
                                     (A0)(A1)(T),
                                     (scalar_<integer_ <A0> )
                                     (scalar_<integer_ <A1> )
                                     (target_<scalar_<floating_ <T> > )
                                     )
  {
    typedef typename T::value_type value_type; 
    typedef table<value_type> result_type; 
    inline result_type operator()(const A0 & n, const A1& k,
                                  const T &) const
      { 
        
        n--;
        result_type  t;
        result_type x = colon(value_type(n), Mone<value_type>(), One<value_type>());
        result_type z = colon(One<value_type>(), One<value_type>(), value_type(n));
        if (k == 0)
          {
            t = diag(x, -1) + diag(z, 1);
          }
        else
          {
            result_type y= sqrt(mul(x,z));
            t = diag(y, -1) + diag(y, 1);
          }; 
        return t;
      };  
  }
}

#endif
