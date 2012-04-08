/*******************************************************************************
 *         Copyright 2003-2012 LASME UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_FUNCTIONS_DRAMADAH_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_FUNCTIONS_DRAMADAH_HPP_INCLUDED
#include <nt2/include/functions/eye.hpp>
#include <nt2/include/functions/zeros.hpp>
#include <nt2/include/functions/colon.hpp>
#include <nt2/include/functions/toeplitz.hpp>

namespace nt2
{
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::dramadah_, tag::cpu_,
                                     (A0)(T), 
                                     (scalar_<integer_<A0> > )
                                     (target_<scalar_<floating_<T> > >)
                                     )
  { 
    typedef typename T::value_type value_type; 
    typedef table< value_type > result_type; 
    inline result_type operator()(const A0 & n, const T&) const
      {
        return dramadah(n, Zero<A0>(), T()); 
      }
  };
  
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::dramadah_, tag::cpu_,
                                     (A0)(A1)(T),
                                     (scalar_<integer_ <A0> )
                                     (scalar_<integer_ <A1> )
                                     (target_<scalar_<floating_<T> > >)
                                     )
  {
    typedef T::value_type value_type; 
    typedef table<value_type, matlab_index_> result_type; 
    inline result_type operator()(const A0 & n, const A1& k,
                                  const T & ) const
      { 
        if(k == 1)
          {
            // Toeplitz
            result_type c = ones(n,1,T());
            for(int i=2; i <= n; i+= 4)
              {
                size_t m = nt2::min(1,n-i);
                c(Range(i, i+m)) = Zero<value_type>();
              }
            result_type r = zeros(n,1,T());
            r(1) = r(2) = r(4) = One<value_type>(); 
            if (n < 4) r = r(Range (1, n));
            return toeplitz (c,r);
          } 
        else if (k == 2)
          { // Upper triangular and Toeplitz
            result_type c = zeros(n,1,T());
            c(1) = One<value_type>();
            result_type r = ones(n,1,T());
            r(Range (3, 2, n)) = Zero<value_type>();
            return toeplitz (c,r);
          } 
        else 
          { //   Lower Hessenberg.
            result_type  c = ones(n,1,T());
            c(Range (2, 2, n)) = T(0);
            result_type  d = zeros(n,1,T());
            d(1) = d(2) = One<value_type>(); 
            return toeplitz(c, d);
          }
      };  
  }
}

#endif
