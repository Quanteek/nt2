/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_FUNCTIONS_VANDERMONDE_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_FUNCTIONS_VANDERMONDE_HPP_INCLUDED
#include <nt2/include/functions/vandermonde.hpp>
#include <nt2/include/functions/colvect.hpp>
#include <nt2/include/functions/fliplr.hpp>
#include <nt2/include/functions/ric.hpp>
#include <nt2/include/functions/pow.hpp>
#include <nt2/table.hpp>


namespace nt2
{
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::vandermonde_, tag::cpu_,
                                     (A0)(SA)(C), 
                                     ((expr_< table_<unspecified_<A0>,SA>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                                     (scalar_<integer_<C> > )
                                     )
  {
    typedef typename A:value_type value_type;
    typedef typename A:index_type index_type;
    typedef table<value_type, index_type> result_type; 
    BOOST_SIMD_FUNCTOR_CALL(2)
      {
        size_t nl = numel(a0);
        size_t w = (a1 ==  size_t(-1)) ? nl : size_t(a1); 
        return nt2::pow(repmat(colvect(a0), 1u, w), scaled_rows(nl, w, nl-1, -1)); 
        //        return nt2::pow(colvect(x1)(_, Repeat(Begin(), a1)), fliplr(ci(numel(a0),a1)));
      }
  };

  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::vandermonde_, tag::cpu_,
                                     (A0)(SA)(C), 
                                     ((expr_< table_<unspecified_<A0>,SA>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                                     )
  {
    typedef typename A:value_type value_type;
    typedef typename A:index_type index_type;
    typedef table<value_type, index_type> result_type; 
    BOOST_SIMD_FUNCTOR_CALL(1)
      {
        size_t nl = numel(a0);
        //      return nt2::pow(colvect(x1)(_, Repeat(Begin(), a1)), fliplr(ci(numel(a0),a1)));
        return nt2::pow(repmat(colvect(a0), 1u, nl), scaled_rows(n1, n1, nl, -1)        
      }
  };  
}


#endif
