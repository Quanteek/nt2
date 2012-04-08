/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_FUNCTIONS_HILBERT_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_FUNCTIONS_HILBERT_HPP_INCLUDED



namespace nt2
{
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::hilbert_, tag::cpu_,
                                     (A0)(T), 
                                     (scalar_<integer_<A1> > )
                                     (target_<floating<T> >)
                                     )
  {
    typedef table< double > result_type; 
    BOOST_SIMD_FUNCTOR_CALL(2)
      {
        typedef typename A1::value_type value_type; 
        return rec(rif(n, T())+ci(n, T()));
        // return bsxfun(tag::plus_, cif(1, n), ri(n, 1));  //Is it better ?
      }
  };
  
}


#endif

// /////////////////////////////////////////////////////////////////////////////
// End of hilbert.hpp
// /////////////////////////////////////////////////////////////////////////////
