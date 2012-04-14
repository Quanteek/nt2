/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_SOLVERS_INTERNAL_SOLVE_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_SOLVERS_INTERNAL_SOLVE_HPP_INCLUDED
#include <nt2/include/functions/solve_lu_ip.hpp>
#include <nt2/include/functions/solve_qr_ip.hpp>
#include <nt2/include/functions/solve_svd_ip.hpp>
#include <nt2/include/functions/solve_chol_ip.hpp>
#include <nt2/include/functions/solve_tr_ip.hpp>
#include <nt2/include/functions/pqr.hpp>
#include <nt2/include/functions/istriangular.hpp>
#include <nt2/include/functions/istriu.hpp>
#include <nt2/include/functions/istril.hpp>
#include <nt2/include/functions/ishermitian.hpp>
#include <iostream>

////////////////////////////////////////////////////////////////////////
// At the exception of ml and mrdivide, 
// all routines listed are called as solve_xxx(a, b) [ solve_xxx(a, b, info)]
// or solve_xxx(a, b, allowdestroy()) [ solve_xxx(a, b, info, allowdestroy())]
// xxx being lu, qr, svd, tril, triu, chol
//
// a is a matrix
// b is a matrix of column right-hand members
// to solve a*x = b and returning x a matrix of solutions
// if allowdestroy() is present the inputs a and b can be used as
// working storage for computations and thus the contained datas
// can be destroyed,  but the counterpart is that somme of the datas
// need not to be copied to insure the preservation.
////////////////////////////////////////////////////////////////////////

namespace nt2
{

  struct allowdestroy {};     
  struct methods
  {
    enum vals { fail    =   0,
                chol    =   1,
                lu      =   2,
                pqr     =   3,
                qr      =   4,
                svd     =   5,
                triu    =  50,
                tril    =  51,
                triperm =  52,
                hess    =  53, 
                spdiag  = 100,
                spsqbd  = 101
    };
    
    methods(vals v = fail):val(v){}
    ~methods(){}
    methods & operator = (const vals & v)
    {
      val = v;
      return *this;
    }
    bool operator == (const methods & m)
    {
      return m.val == val;
    }
    vals val;     
  }; 
  
  inline std::ostream & operator << (std::ostream & c, const methods & m){
    switch (m.val)
      {
      case methods::fail     : c << "fail";   break; 
      case methods::chol     : c << "chol";   break;
      case methods::lu       : c << "lu";     break;  
      case methods::pqr      : c << "pqr";    break; 
      case methods::qr       : c << "qr";     break;
      case methods::svd      : c << "svd";    break;
      case methods::triu     : c << "triu";   break;
      case methods::tril     : c << "tril";   break;
      case methods::triperm  : c << "triperm";break;
      case methods::hess     : c << "hess";   break;
      case methods::spdiag   : c << "spdiag"; break;
      case methods::spsqbd   : c << "spsqbd"; break;
      }
    return c; 
  }

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
    
    static inline tab_t solve_lu(const A &a, const B &b, la_int & info)
    {
      tab_t aa(a);
      tab_t bb(b);
      info = typename nt2::make_functor<tag::solve_lu_ip_, tab_t>::type()(aa, bb).get_info();
      return bb; 
    }
    static inline tab_t solve_lu_ip(A &a, B &bx, la_int & info)
    {
      info = typename nt2::make_functor<tag::solve_lu_ip_, A>::type()(a, bx).get_info();
      return bx; 
    }
    static inline tab_t solve_qr(const A &a, const B &b, la_int & info)
    {
      tab_t aa(a);
      tab_t x;
      info = typename nt2::make_functor<tag::solve_qr_ip_, tab_t>::type()(aa, x, b).get_info();
      return x; 
    }
    static inline tab_t solve_qr_ip(A &a, B &bx, la_int & info)
    {
      info = typename nt2::make_functor<tag::solve_qr_ip_, A>::type()(a, bx, bx).get_info();
      return bx; 
    }
    static inline tab_t solve_pqr_ip(A &a, B &b, la_int & info)
    {
      nt2::pqr_return<A> f = pqr(a);
      return f.solve(b); ; 
    }
     static inline tab_t solve_pqr(const A &a, const B &b, la_int & info)
     {
      tab_t aa(a);
      tab_t bb(b);
      nt2::pqr_return<tab_t> f = pqr(aa);
      info = f.get_info();
      return f.solve(bb); 
     }
    static inline tab_t solve_svd(const A &a, const B &b, la_int & info)
    {
      tab_t aa(a);
      tab_t x;
      info = typename nt2::make_functor<tag::solve_svd_ip_, tab_t>::type()(aa, x, b).get_info();
      return x; 
    }
    static inline tab_t solve_svd_ip(A &a, B &bx, la_int & info)
    {
      info = typename nt2::make_functor<tag::solve_svd_ip_, A>::type()(a, bx, bx).get_info();
      return bx; 
    }

    static inline tab_t solve_chol(const A &a, const B &b, la_int & info)
    {
      tab_t aa(a);
      tab_t bb(b);
//       nt2::solve_chol_ip_return<tab_t, tab_t> f = solve_chol_ip(aa, bb, 'l');
//       info = f.get_info(); 
      info = typename nt2::make_functor<tag::solve_chol_ip_, tab_t>::type()(aa, bb, 'l').get_info();
      return bb; 
    }
    static inline tab_t solve_chol_ip(A &a, B &bx, la_int & info)
    {
      info = typename nt2::make_functor<tag::solve_chol_ip_, tab_t>::type()(a, bx, 'l').get_info();
      return bx; 
    }


    static inline tab_t solve_tril(const A &a, const B &b, la_int & info)
    {
      tab_t bb(b);
      info = typename nt2::make_functor<tag::solve_tr_ip_, tab_t>::type()(a, bb, 'l', 'n', 'n').get_info();
      return bb; 
    }
    static inline tab_t solve_tril_ip(A &a, B &bx, la_int & info)
    {
      info = typename nt2::make_functor<tag::solve_tr_ip_, A>::type()(a, bx, 'l', 'n', 'n').get_info();
      return bx; 
    }
    static inline tab_t solve_triu(const A &a, const B &b, la_int & info)
    {
      tab_t bb(b);
      info = typename nt2::make_functor<tag::solve_tr_ip_, tab_t>::type()(a, bb, 'u', 'n', 'n').get_info();
      return bb; 
    }
    static inline tab_t solve_triu_ip(A &a, B &bx, la_int & info)
    {
      info = typename nt2::make_functor<tag::solve_tr_ip_, A>::type()(a, bx, 'u', 'n', 'n').get_info();
      return bx; 
    }
    static inline tab_t mldivide(const A &a, const B &b, la_int & info, methods& meth)
    {
      tab_t x;
      meth =  methods::fail;
      info = 0; 
      if (is_triu(a))
        {
          std::cout <<  "triu" << std::endl; 
          x = solve_triu(a,b,info);
          std::cout << "info " << info << std::endl; 
          if (info == 0) meth = methods::triu; 
        }
      if (meth == methods::fail && is_tril(a))
        {
          std::cout <<  "tril" << std::endl; 
          x = solve_tril(a,b,info);
          std::cout << "info " << info << std::endl; 
          if (info == 0) meth = methods::tril;
        }
      if (meth == methods::fail && is_hermitian(a))
        {
          std::cout <<  "chol" << std::endl; 
          x = solve_chol(a,b,info);
          std::cout << "info " << info << std::endl; 
          if(info == 0) meth = methods::chol;
        }
      if (meth == methods::fail && is_square(a))
        {
          std::cout <<  "lu" << std::endl; 
          x = solve_lu(a, b,info); 
          std::cout << "info " << info << std::endl; 
          if(info == 0) meth = methods::lu; 
        }
      if (meth == methods::fail)
        { 
          std::cout <<  "pqr" << std::endl; 
          x = solve_pqr(a,b,info);
          std::cout << "info " << info << std::endl; 
          meth = methods::pqr; 
        }
      std::cout << "method " << meth << "  info =  " << info << std::endl; 
      return x; 
    }    
//     static inline tab_t mrdivide(const B &b, const A &a, la_int & info, methods& meth)
//     {
//       return trans(mldivide(trans(a), trans(b), info, meth));              
//     }
  }; 
  
#define NT2_SOLVER(S)                                           \
  template < class A, class B> typename solve<A, B>::tab_t      \
  BOOST_PP_CAT(solve_, S)(const A &a, const B &b, la_int& info) \
  {                                                             \
    return solve<A, B>::BOOST_PP_CAT(solve_,S)(a,b,info);       \
  }                                                             \
  template < class A, class B> typename solve<A, B>::tab_t      \
  BOOST_PP_CAT(solve_, S)(A &a, B &b, la_int& info, const allowdestroy &) \
  {                                                                     \
    return solve<A, B>::BOOST_PP_CAT(BOOST_PP_CAT(solve_,S),_ip)(a,b,info); \
  }                                                                     \
  template < class A, class B> typename solve<A, B>::tab_t            \
  BOOST_PP_CAT(solve_, S)(const A &a, const B &b)                       \
  {                                                                     \
    la_int  info;                                                       \
    return solve<A, B>::BOOST_PP_CAT(solve_,S)(a,b,info);               \
  }                                                                     \
  template < class A, class B> typename solve<A, B>::tab_t            \
  BOOST_PP_CAT(solve_, S)(A &a, B &b, const allowdestroy &)             \
  {                                                                     \
    la_int  info;                                                       \
    return solve<A, B>::BOOST_PP_CAT(BOOST_PP_CAT(solve_,S),_ip)(a,b,info); \
  }                                                                     \
    
  NT2_SOLVER(pqr)
  NT2_SOLVER(lu)
  NT2_SOLVER(qr)
  NT2_SOLVER(svd)
  NT2_SOLVER(chol)
  NT2_SOLVER(tril)
  NT2_SOLVER(triu)
#undef NT2_SOLVER
    
  template < class A, class B> typename solve<A, B>::tab_t      
  mldivide(const A &a, const B &b, methods& meth) 
  {
    la_int info; 
    typename solve<A, B>::tab_t x = solve<A, B>::mldivide(a, b, info, meth);
    BOOST_ASSERT_MSG(info == 0, "mldivide failed"); 
    return x; 
  }                                                             
  
}


#endif
