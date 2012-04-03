/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_FUNCTIONS_MAGIC_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_FUNCTIONS_MAGIC_HPP_INCLUDED



namespace nt2
{
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::magic_, tag::cpu_,
                                     (A0), 
                                     (scalar_<integer_<A0> > )
                                     )
  {
    typedef table< ptrdiff_t > result_type; 
    inline result_type operator()()
      {
        //  The matrix is 8-by-8 with integer elements.
        //  It has:
        //      * A double eigenvalue.
        //      * Three nearly equal eigenvalues.
        //      * Dominant eigenvalues of opposite sign.
        //      * A zero eigenvalue.
        //      * A small, nonzero eigenvalue.
        static ptrdiff_t r[] = {
          611,  196, -192,  407, -8,  -52,  -49,   29, 
          196,  899,  113, -192,  -71,  -43,   -8,  -44, 
          -192,  113,  899,  196, 61,   49,    8,   52, 
          407, -192,  196,  611,    8,   44,   59,  -23, 
          -8,  -71,   61,    8,  411, -599,  208,  208, 
          -52,  -43,   49,   44, -599,  411,  208,  208, 
          -49,  -8,    8,   59,  208,  208,   99, -911, 
          29,  -44,   52,  -23,   208,  208,  -911,   99};
        return table < ptrdiff_t > (&r[0], of_size_<8, 8>); 
    }
  };
  
}


#endif

// /////////////////////////////////////////////////////////////////////////////
// End of magic.hpp
// /////////////////////////////////////////////////////////////////////////////
