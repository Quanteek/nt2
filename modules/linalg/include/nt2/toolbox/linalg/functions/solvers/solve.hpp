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
    
  }; 

  template < class A, class B> typename solve<A, B>::tab_t
  solve_lu(const A &a, const B &b)
  {
    return solve<A, B>::solve_lu(a,b); 
  }
  template < class A, class B> typename solve<A, B>::tab_t
  solve_lu(A &a, B &b, const allowdestroy &)
  {
    return solve<A, B>::solve_lu_ip(a,b); 
  }
  template < class A, class B> typename solve<A, B>::tab_t
  solve_qr(const A &a, const B &b)
  {
    return solve<A, B>::solve_qr(a,b); 
  }
  template < class A, class B> typename solve<A, B>::tab_t
  solve_qr(A &a, B &b, const allowdestroy &)
  {
    return solve<A, B>::solve_qr_ip(a,b); 
  }
  template < class A, class B> typename solve<A, B>::tab_t
  solve_svd(const A &a, const B &b)
  {
    return solve<A, B>::solve_svd(a,b); 
  }
  template < class A, class B> typename solve<A, B>::tab_t
  solve_svd(A &a, B &b, const allowdestroy &)
  {
    return solve<A, B>::solve_svd_ip(a,b); 
  }    
}


#endif

// /////////////////////////////////////////////////////////////////////////////
// End of solve_lu.hpp
// /////////////////////////////////////////////////////////////////////////////
