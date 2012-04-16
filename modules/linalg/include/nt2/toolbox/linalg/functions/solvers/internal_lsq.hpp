/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_SOLVERS_INTERNAL_LSQ_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_SOLVERS_INTERNAL_LSQ_HPP_INCLUDED
#include <nt2/include/functions/lsq_lse_ip.hpp>
#include <iostream>

////////////////////////////////////////////////////////////////////////
// all routines listed are called as lsq_xxx(a, b) [ lsq_xxx(a, b, info)]
// or lsq_xxx(a, b, allowdestroy()) [ lsq_xxx(a, b, info, allowdestroy())]
// xxx being lu, qr, svd, tril, triu, chol
//
// a is a matrix
// b is a matrix of column right-hand members
// to lsq a*x = b and returning x a matrix of solutions
// if allowdestroy() is present the inputs a and b can be used as
// working storage for computations and thus the contained datas
// can be destroyed,  but the counterpart is that somme of the datas
// need not to be copied to insure the preservation.
////////////////////////////////////////////////////////////////////////

namespace nt2
{

  template < class A, class B, class C, class D, class X> struct lsq {
    
    typedef typename A::value_type                   type_t;
    typedef typename A::index_type                  index_t; 
    typedef typename meta::as_real<type_t>::type    btype_t; 
    typedef nt2::table<type_t,nt2::matlab_index_>    ftab_t;
    typedef nt2::table<btype_t,nt2::matlab_index_>  fbtab_t;
    typedef nt2::table<la_int,nt2::matlab_index_>   fitab_t;
    typedef nt2::table<type_t,index_t>                tab_t;
    typedef nt2::table<btype_t,index_t>              btab_t;
    typedef nt2::table<la_int,index_t>               itab_t;
    
    static inline tab_t lsq_lse(const A &a, const B &b,const C &c, const D &d, X& x, la_int & info)
    {
      tab_t aa(a);
      tab_t bb(b);
      tab_t aa(c);
      tab_t bb(d);
      info = typename nt2::make_functor<tag::lsq_lse_ip_, tab_t>::type()(aa, bb, cc, dd, x).get_info();
      return x; 
    }
    static inline tab_t lsq_lse_ip(A &a, B &b, C &c, D &d, X& x, la_int & info)
    {
      info = typename nt2::make_functor<tag::lsq_lu_ip_, A>::type()(a, b, c, d, x).get_info();
      return x; 
    }
  }; 
  
#define NT2_LSQ(S)                                                      \
  template < class A, class B> typename lsq<A, B>::tab_t                \
  BOOST_PP_CAT(lsq_, S)(const A &a, const B &b, la_int& info)           \
  {                                                                     \
    return lsq<A, B>::BOOST_PP_CAT(lsq_,S)(a,b,info);                   \
  }                                                                     \
  template < class A, class B> typename lsq<A, B>::tab_t                \
  BOOST_PP_CAT(lsq_, S)(A &a, B &b, la_int& info, const allowdestroy &) \
  {                                                                     \
    return lsq<A, B>::BOOST_PP_CAT(BOOST_PP_CAT(lsq_,S),_ip)(a,b,info); \
  }                                                                     \
  template < class A, class B> typename lsq<A, B>::tab_t                \
  BOOST_PP_CAT(lsq_, S)(const A &a, const B &b)                         \
  {                                                                     \
    la_int  info;                                                       \
    return lsq<A, B>::BOOST_PP_CAT(lsq_,S)(a,b,info);                   \
  }                                                                     \
  template < class A, class B> typename lsq<A, B>::tab_t                \
  BOOST_PP_CAT(lsq_, S)(A &a, B &b, const allowdestroy &)               \
  {                                                                     \
    la_int  info;                                                       \
    return lsq<A, B>::BOOST_PP_CAT(BOOST_PP_CAT(lsq_,S),_ip)(a,b,info); \
  }                                                                     \
    
  NT2_LSQ(lse)

#undef NT2_LSQ
  
}


#endif
