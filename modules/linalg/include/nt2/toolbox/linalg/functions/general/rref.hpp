/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_GENERAL_RREF_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_GENERAL_RREF_HPP_INCLUDED

namespace nt2
{
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::chol_, tag::cpu_, 
                              (A)(S0)(A1),
                              ((expr_< table_<unspecified_<A>,S0>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                              (scalar_<floating_<A1> > )
                              )
  {
    typedef typename A::value_type value_type;
    typedef typename meta::as_integer<value_type, unsigned> itype_t; 
    typedef nt2::table<value_type> tab_t;
    typedef nt2::table<itype_t>    itab_t;
    typedef boost::mpl_vector < tab_t, itab_t > result_type; 
    BOOST_DISPATCH_FORCE_INLINE result_type operator()(const A& aa, const A1& tol) const
    {
      tab_t a = aa; 
      itype_t i = 1, j = 1; 
      itype_t n = size(a, 1);
      itype_t m = size(a, 2);
      itype_t p, k; 
      while(i <= m && j <= n)
        {
          tie(p, k) =  nt2::max(nt2::abs(a(_(i, m),j)));
          k = k+i-1;
          if (p <= tol)
            {
            // the column is negligible, zero it out.
              a(_(i, m),j) = Zero<value_type>();
              ++j;
            }
          else
            {
              // remember column index
              jb = cath(jb, j);
              // swap i-th and k-th rows.
              a(cath(i, k),_(j, n)) = a(cath(k, i),_(j, n));
              // divide the pivot row by the pivot element.
              a(i,cath(j, n)) = a(i,cath(j, n))/a(i,j);
              // subtract multiples of the pivot row from all the other rows.
              for (itype_t k = 1; k <= m; ++k)//[1:i-1 i+1:m]
                {
                  if (k!=i) {
                    a(k,_(j, n)) = a(k,_(j, n)) - a(k,j)*a(i,_(j, n));
                  }
                }
              ++i; ++j;
            }
        }
      return result_type(a, jb); 
    }
  };

  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::chol_, tag::cpu_, 
                              (A)(S0)(A1),
                              ((expr_< table_<unspecified_<A>,S0>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                              (scalar_<floating_<A1> > )
                              )
  {
    typedef typename A::value_type value_type; 
    typedef nt2::table<value_type> result_type; 
    BOOST_DISPATCH_FORCE_INLINE result_type operator()(const A& a) const
    {
      value_type tol = nt2::max(size(a))*nt2::eps(nt2::norm(A,'inf'));
      return rref(a, tol); 
    }
  };  
}


#endif

// /////////////////////////////////////////////////////////////////////////////
// End of rref.hpp
// /////////////////////////////////////////////////////////////////////////////
