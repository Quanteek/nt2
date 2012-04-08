/*******************************************************************************
 *         Copyright 2003-2012 LASME UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_FUNCTIONS_CYCOL_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_FUNCTIONS_CYCOL_HPP_INCLUDED
#include <nt2/include/functions/eye.hpp>
#include <nt2/include/functions/zeros.hpp>
#include <nt2/include/functions/colon.hpp>
#include <nt2/include/functions/toeplitz.hpp>
#include <nt2/include/functions/pow.hpp>

namespace nt2
{
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::cycol_, tag::cpu_,
                                     (A0)(A1)(T), 
                                     (scalar_<integer_<A0> > )
                                     (scalar_<integer_<A1> > )
                                     (target_<scalar_<floating_<T> > >)
                                     )
  { 
    typedef typename T::value_type value_type; 
    typedef table< value_type > result_type; 
    inline result_type operator()(const A0 & m, const A0 & n, const T&) const
      {
        return cycol(m, n, 0, T()); 
      }
  };
  
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::cycol_, tag::cpu_,
                                     (A0)(A1)(A2)(T),
                                     (scalar_<integer_ <A0> )
                                     (scalar_<integer_ <A1> )
                                     (scalar_<integer_ <A2> )
                                     (target_<scalar_<floating_ <T> > )
                                     )
  {
    typedef typename T::value_type value_type; 
    typedef table<value_type> result_type; 
    inline result_type operator()(const A0 & m, const A1& n,const A2& k0,
                                  const T &) const
      { 
        if(k0 == 0) k = size_t(nt2::max(idivround(n, Four<A1>()),One<A2>()));
      result_type a = randn(m, k); 
      result_type c(of_size(m, 0)),  c1; 
      size_t deb = 0;
      for (size_t i=1;  i <= nt2::idivceil(n, k);  i++){
        c1 = cath(c, a(_,AllTo(k-1)));
        c =  c1; 
        deb+= k; 
      }
     a = c(_, AllTo(n-1));
      return a;
      };  
  }
}

#endif
