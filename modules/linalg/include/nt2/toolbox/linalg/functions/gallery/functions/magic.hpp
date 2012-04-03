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
    BOOST_SIMD_FUNCTOR_CALL(1)
      {
        // %MAGIC  Magic square.
        // %   MAGIC(a0) is an a0-by-a0 table constructed from the integers
        // %   1 through a0^2 with equal row, column, and diagonal sums.
        // %   Produces valid magic squares for all a0 > 0 except a0 = 2.
    
      if (a0%2 == 1)   // % Odd order.
        {
          table < double > I, J; 
          meshgrid(nt2::column(1, a0), nt2::column(1, a0), J, I);
          table < double > A = nt2::mod(I+J-(a0+3.0)/2.0,a0);
          table < double > B = nt2::mod(I+2.0*J-2.0,a0);
          return nt2::convert < ptrdiff_t > (double(a0)*A + B + 1.0);
        }
      else if (a0%4 == 0)   // % Doubly even order.
        {
          table < ptrdiff_t > I, J; 
          meshgrid(nt2::column(1, a0), nt2::column(1, a0), J, I);
          table < bool > K = fix(mod(I,4)/2) == fix(mod(J,4)/2);
          table < double > M = trans(reshape(nt2::column(1, sqr(a0)),a0,a0));
          M(K) = sqr(a0) + 1 - M(K);
          return nt2::convert < ptrdiff_t > (M); 
        }
      else  // % Singly nt2n order.
        {
          size_t  p = a0/2;
          table < double > M = magic(p);
          M = nt2::catv( nt2::cat(M)(M+2*sqr(p))(),
                         nt2::cat(M+3*sqr(p))(M+sqr(p))());
          if (a0 == 2) return  nt2::convert < ptrdiff_t > (M); 
          table < ptrdiff_t > i = colvect(nt2::column(1, p));
          size_t  k = (a0-2)/4;
          table < ptrdiff_t > j = nt2::cath(nt2::column(1, k), nt2::column((a0-k+2), a0));
          M(nt2::catv(i, i+p),j) = M(nt2::catv(i+p, i),j);
          i = nt2::cons(k+1);
          j = nt2::cat(1)(i)();
          M(nt2::catv(i, i+p),j) = M(nt2::catv(i+p, i),j);
          return  nt2::convert < ptrdiff_t > (M);
        }
      return zeros(a0); 
    }
  };
  
}


#endif

// /////////////////////////////////////////////////////////////////////////////
// End of magic.hpp
// /////////////////////////////////////////////////////////////////////////////
