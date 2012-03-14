/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_SOLVERS_SOLVE_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_SOLVERS_SOLVE_HPP_INCLUDED
#include <nt2/include/functions/solve_lu_ip.hpp>
#include <nt2/include/functions/solve_qr_ip.hpp>
#include <nt2/include/functions/solve_svd_ip.hpp>
#include <nt2/include/functions/solve_chol_ip.hpp>
#include <nt2/include/functions/solve_tr_ip.hpp>

////////////////////////////////////////////////////////////////////////
// all routines listed are called as solve_xxx(a, b)
// or solve_xxx(a, b, allowdestroy())
// xxx being lu, qr, svd, tril, triu, chol
//
// a is a matrix
// b is a matix of second members
// to solve a*x = b and returning x
// if allowdestroy() is present the input a and b can be used as
// working storage for computation and thus the contained datas
// can be destroyed.
////////////////////////////////////////////////////////////////////////

namespace nt2
{

  struct allowdestroy {};     

  template < class A, class B> struct solve {
    
    typedef typename A::value_type                   type_t;
    typedef typename A::index_type                  index_t; 
    typedef typename meta::as_real<type_t>::type    btype_t; 
    typedef nt2::table<type_t,nt2::matlab_index_>    ftab_t;
    typedef nt2::table<btype_t,nt2::matlab_index_>  fbtab_t;
    typedef nt2::table<la_int,nt2::matlab_index_>   fitab_t;
    typedef nt2::table<type_t,index_t>                tab_t;
    typedef nt2::table<btype_t,index_t>              btab_t;
    typedef nt2::table<la_int,index_t>               itab_t;
    
    static inline tab_t solve_lu(const A &a, const B &b)
    {
      tab_t aa(a);
      tab_t bb(b);
      typename nt2::make_functor<tag::solve_lu_ip_, tab_t>::type()(aa, bb);
      return bb; 
    }
    static inline tab_t solve_lu_ip(A &a, B &bx)
    {
      typename nt2::make_functor<tag::solve_lu_ip_, A>::type()(a, bx);
      return bx; 
    }
    static inline tab_t solve_qr(const A &a, const B &b)
    {
      tab_t aa(a);
      tab_t x;
      typename nt2::make_functor<tag::solve_qr_ip_, tab_t>::type()(aa, x, b);
      return x; 
    }
    static inline tab_t solve_qr_ip(A &a, B &bx)
    {
      typename nt2::make_functor<tag::solve_qr_ip_, A>::type()(a, bx, bx);
      return bx; 
    }
    static inline tab_t solve_svd(const A &a, const B &b)
    {
      tab_t aa(a);
      tab_t x;
      typename nt2::make_functor<tag::solve_svd_ip_, tab_t>::type()(aa, x, b);
      return x; 
    }
    static inline tab_t solve_svd_ip(A &a, B &bx)
    {
      typename nt2::make_functor<tag::solve_svd_ip_, A>::type()(a, bx, bx);
      return bx; 
    }
    static inline tab_t solve_chol(const A &a, const B &b)
    {
      tab_t aa(a);
      tab_t bb(b);
      typename nt2::make_functor<tag::solve_chol_ip_, tab_t>::type()(aa, bb);
      return bb; 
    }
    static inline tab_t solve_chol_ip(A &a, B &bx)
    {
      typename nt2::make_functor<tag::solve_chol_ip_, A>::type()(a, bx);
      return bx; 
    }
    static inline tab_t solve_tril(const A &a, const B &b)
    {
      tab_t bb(b);
      typename nt2::make_functor<tag::solve_tr_ip_, tab_t>::type()(a, bb, 'l', 'n', 'n');
      return bb; 
    }
    static inline tab_t solve_tril_ip(A &a, B &bx)
    {
      typename nt2::make_functor<tag::solve_tr_ip_, A>::type()(a, bx, 'l', 'n', 'n');
      return bx; 
    }
     static inline tab_t solve_triu(const A &a, const B &b)
    {
      tab_t bb(b);
      typename nt2::make_functor<tag::solve_tr_ip_, tab_t>::type()(a, bb, 'u', 'n', 'n');
      return bb; 
    }
    static inline tab_t solve_triu_ip(A &a, B &bx)
    {
      typename nt2::make_functor<tag::solve_tr_ip_, A>::type()(a, bx, 'u', 'n', 'n');
      return bx; 
    }    
  }; 

#define NT2_SOLVER(S)                                           \
  template < class A, class B> typename solve<A, B>::tab_t      \
  BOOST_PP_CAT(solve_, S)(const A &a, const B &b)               \
  {                                                             \
    return solve<A, B>::BOOST_PP_CAT(solve_,S)(a,b);            \
  }                                                             \
  template < class A, class B> typename solve<A, B>::tab_t      \
  BOOST_PP_CAT(solve_, S)(A &a, B &b, const allowdestroy &)     \
  {                                                             \
    return solve<A, B>::BOOST_PP_CAT(BOOST_PP_CAT(solve_,S),_ip)(a,b); \
  }                                                             \

  NT2_SOLVER(lu)
  NT2_SOLVER(qr)
  NT2_SOLVER(svd)
  NT2_SOLVER(chol)
  NT2_SOLVER(tril)
  NT2_SOLVER(triu)  
 #undef NT2_SOLVER 

}


#endif

// /////////////////////////////////////////////////////////////////////////////
// End of solve_lu.hpp
// /////////////////////////////////////////////////////////////////////////////
