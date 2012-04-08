/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_FUNCTIONS_CHEBVAND_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_FUNCTIONS_CHEBVAND_HPP_INCLUDED
#include <nt2/include/functions/linspace.hpp>
#include <nt2/include/functions/ones.hpp>

namespace nt2
{
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::chebvand_, tag::cpu_,
                                     (A0)(T), 
                                     (scalar_<integer_<A1> > )
                                     (target_<floating<T> >)
                                     )
  { 
    typedef typename T:::value_type value_type; 
    typedef table< value_type > result_type; 
    BOOST_SIMD_FUNCTOR_CALL(2)
      {
        typedef typename A1::value_type                  value_type; 
        typedef typename table<value_type, matlab_index_>   table_t; 
        table_t p = linspace(value_type(0),value_type(1),n);
        table_t c = ones(n,n,meta::as_<value_type>());
        if (n == 1) return c;
        c(1, _) = p;
        //      Use Chebyshev polynomial recurrence.
        for (size_t i = 2; i < n; i++){
          c(i, _) = 2*mul(p, c(i-1, _)-c(i-2, _));
        }
        return c;
        
      }
  };
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::chebvand_, tag::cpu_,
                                     (A0)(S0)(A1),
                                     ((expr_< table_<unspecified_<A0>,S0>,nt2::tag::terminal_,boost::mpl::long_<0> >)),
                                     (scalar_<integer_ < A1> > )
                                     )
  { 
    typedef typename A0:::value_type value_type; 
    typedef table< value_type > result_type; 
    inline result_type operator()(const A0&p, const A1&mm))
      {
        typedef typename A1::value_type   value_type; 
        typedef typename table<value_type, matlab_index_>   table_t; 
        size_t n = numel(p); 
        size_t m = mm ? mm : n;
        matrix < T > c = ones(m,n);
        if (m == 1) return c; 
        c(1, _) = p; 
        //      Use Chebyshev polynomial recurrence.
        for (size_t i = 2; i < m; i++){
          c(i, _)) = 2*mul(p,c(i-1, _))-c(i-2, _);
        }
        return c; 
      }
  };  
}


#endif

// /////////////////////////////////////////////////////////////////////////////
// End of chebvand.hpp
// /////////////////////////////////////////////////////////////////////////////
