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
                                     (A0)(T), 
                                     (scalar_<integer_<A0> > )
                                     ((target_ < scalar_<floating_<T> > > ))
                                     )
  {
    typedef T::value_type value_type; 
    typedef table< value_type > result_type; 
    BOOST_SIMD_FUNCTOR_CALL(2)
      {
        BOOST_ASSERT_MSG(n!= 2, "n is an integer no equal to 2"); 
        if (a0%2 == 1)   // % Odd order.
          {
            result_type I = rif(1, a0, T());
            result_type J = cif(1, a0, T()); 
            result_type A = nt2::mod(I+J-(a0+3.0)/2.0,a0);
            result_type B = nt2::mod(I+2.0*J-2.0,a0);
            return value_type(a0)*A + B + One<value_type>();
          }
        else if (a0%4 == 0)   // % Doubly even order.
          {
            table < ptrdiff_t > I, J; 
            result_type I = rif(1, a0, meta::as_<value_type>());
            result_type J = cif(1, a0, meta::as_<value_type>()); 
            table < logical<value_type> > K = nt2::fix(nt2::mod(I,4)/2) == nt2::fix(nt2::mod(J,4)/2);
            result_type M = trans(reshape(nt2::column(1, sqr(a0), T()),a0,a0));
            M(K) = oneplus(sqr(value_type(a0))) - M(K);
            return nt2::convert < ptrdiff_t > (M); 
          }
        else  // % Singly nt2n order.
          {
            size_t  p = a0/2;
            result_type M = magic(p);
            M = nt2::catv( nt2::cat(M)(M+2*sqr(p))(),
                           nt2::cat(M+3*sqr(p))(M+sqr(p))());
            if (a0 == 2) return  nt2::convert < ptrdiff_t > (M); 
            table < ptrdiff_t > i = colvect(nt2::colon(1, p));
            size_t  k = (a0-2)/4;
            table < ptrdiff_t > j = nt2::cath(nt2::colon(1, k), nt2::colon((a0-k+2), a0));
            M(nt2::catv(i, i+p),j) = M(nt2::catv(i+p, i),j);
            i = nt2::cons(k+1);
            j = nt2::cat(1)(i)();
            M(nt2::catv(i, i+p),j) = M(nt2::catv(i+p, i),j);
            return  M;
          }
        return zeros(a0); 
      }
  };
  
}


#endif

// /////////////////////////////////////////////////////////////////////////////
// End of magic.hpp
// /////////////////////////////////////////////////////////////////////////////
