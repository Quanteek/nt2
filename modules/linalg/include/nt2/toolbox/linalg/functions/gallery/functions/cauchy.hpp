/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_FUNCTIONS_CAUCHY_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_FUNCTIONS_CAUCHY_HPP_INCLUDED
#include <nt2/include/functions/rowvect.hpp>
#include <nt2/include/functions/colvect.hpp>
#include <nt2/include/functions/colon.hpp>
//#include <nt2/include/functions/bsxfun.hpp>



namespace nt2
{
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::cauchy_, tag::cpu_,
                                     (A0)(S0)(A1)(S1),
                                     ((expr_< table_<unspecified_<A0>,S0>,nt2::tag::terminal_,boost::mpl::long_<0> >)), 
                                     ((expr_< table_<unspecified_<A1>,S1>,nt2::tag::terminal_,boost::mpl::long_<0> >)), 
                                     )
  {
    typedef typename A0::value_type value_type; 
    typedef table<value_type, S0> result_type; 
    BOOST_SIMD_FUNCTOR_CALL(2)
      {
        BOOST_ASSERT_MSG(numel(a0) == numel(y), "a0 and y have not the same number of elements"); 
        return rec(bsxfun(nt2::tag::plus_, rowvect(a1), colvect(a0)));
      }
    
  };
  
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::cauchy_, tag::cpu_,
                                     (A0)(A1)(S1),
                                     (scalar_<integer_<A0> >), 
                                     ((expr_< table_<unspecified_<A1>,S1>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                                     )
  {
    typedef typename A0::value_type value_type; 
    typedef table<value_type, S0> result_type; 
    BOOST_SIMD_FUNCTOR_CALL(2)
      {      
        return rec(bsxfun(nt2::tag::plus_, rowvect(_(1, a0)), colvect(a1)));
      }
    
  };
  
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::cauchy_, tag::cpu_,
                                     (A0)(S0)(A1),
                                     ((expr_< table_<unspecified_<A0>,S0>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                                     (scalar_<integer_<A1> >)  
                                     )
  {
    typedef typename A0::value_type value_type; 
    typedef table<value_type, S0> result_type; 
    BOOST_SIMD_FUNCTOR_CALL(2)
      {
        return rec(bsxfun(nt2::tag::plus_, rowvect(a0), colvect(_(1, a1))));
      }
    
  };
  
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::cauchy_, tag::cpu_,
                                     (A0)(A1)(T),
                                     (scalar_<integer_<A0> >)  
                                     (scalar_<integer_<A1> >)
                                     (target_<scalar_<integer_<T> > > )
                                     
                                     )
  {
    typedef typename T::value_type value_type; 
    typedef table<value_type> result_type; 
    BOOST_SIMD_FUNCTOR_CALL(3)
      {
        return rec(bsxfun(nt2::tag::plus_, rowvect(a0), colvect(colon(1, a1))));
      }
    
  };
  
  
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::cauchy_, tag::cpu_,
                                     (A0)(S0),
                                     ((expr_< table_<unspecified_<A0>,S0>,nt2::tag::terminal_,boost::mpl::long_<0> >)), 
                                     )
  {
    typedef typename A0::value_type value_type; 
    typedef table<value_type, S0> result_type; 
    BOOST_SIMD_FUNCTOR_CALL(1)
      {
        return rec(bsxfun(nt2::tag::plus_, rowvect(a0), colvect(a0)));
      }
    
  };
  
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::cauchy_, tag::cpu_,
                                     (A0),
                                     (scalar_<integer_<A0>), 
                                     )
  {
    typedef typename A0::value_type value_type; 
    typedef table<value_type, S0> result_type; 
    BOOST_SIMD_FUNCTOR_CALL(1)
      {
        return rec(bsxfun(nt2::tag::plus_, colon(1, a0), colvect(colon(1, a0))));
      }
    
  };
  
}


#endif
