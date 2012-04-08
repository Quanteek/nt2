/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_FUNCTIONS_FIELDER_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_FUNCTIONS_FIELDER_HPP_INCLUDED



namespace nt2
{
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::fielder_, tag::cpu_,
                                     (A0)(T), 
                                     (scalar_<integer_<A1> > )
                                     (target_<floating<T> >)
                                     )
  {
    typedef typename T::value_type value_type; 
    typedef table< value_type > result_type; 
    BOOST_SIMD_FUNCTOR_CALL(2)
      {
        typedef typename A1::value_type value_type; 
        return fielder(colon(Ones<value_type>(), value_type(c))
      }
  };
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::fielder_, tag::cpu_,
                                     (A0)(SA0), 
                                     (expr_<table_ < unspecified_<A0>, SA0 > > )
                                     )
  {
    typedef typename A0::value_type value_type; 
    typedef table< value_type > result_type; 
    BOOST_SIMD_FUNCTOR_CALL(2)
      {
      result_type d = rowvect(c);//                    % Ensure c is a row vector.
      //result_type a  = nt2::abs(d(Repeat(Begin(), n), _)-colvect(c)(Repeat(Begin(), n), _));
      a =  nt2::abs(bsxfun(tag::minus_, rowvect(c), colvect(c))); 
      return a; 
      }
  };  
}


#endif

// /////////////////////////////////////////////////////////////////////////////
// End of fielder.hpp
// /////////////////////////////////////////////////////////////////////////////
