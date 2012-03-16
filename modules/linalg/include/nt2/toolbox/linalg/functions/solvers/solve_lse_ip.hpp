/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_SOLVERS_SOLVE_LSE_IP_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_SOLVERS_SOLVE_LSE_IP_HPP_INCLUDED

#include <nt2/include/functions/solve_lse_ip.hpp>
#include <nt2/include/functions/leading_size.hpp>
#include <nt2/toolbox/linalg/details/lapack/gelsd.hpp>
#include <nt2/table.hpp>
#include <nt2/include/functions/max.hpp>
#include <nt2/include/functions/min.hpp>
#include <nt2/include/functions/areofsameheight.hpp>
#include <nt2/include/functions/height.hpp>
#include <nt2/include/functions/width.hpp>
//TO DO !!
// #include <nt2/include/functions/diag.hpp>
// #include <nt2/include/functions/expand.hpp>
// #include <nt2/include/functions/fliplr.hpp>
// #include <nt2/include/functions/globalsum.hpp>
// #include <nt2/include/functions/trans.hpp>
// #include <nt2/include/functions/range.hpp>
// #include <nt2/include/functions/first_index.hpp>
// #include <nt2/include/functions/last_height_index.hpp>
#include <iostream>
#include <vector>


//==============================================================================
// lse actual functor forward declaration
//==============================================================================
namespace nt2
{
  template<class A, class X, class B> struct solve_lse_ip_return;
} 

namespace nt2 { namespace ext
{
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::solve_lse_ip_, tag::cpu_, 
                              (A)(SA)(X)(SX)(B)(SB)(C)(SC)(D)(SD),
                              ((expr_< table_<unspecified_<A>,SA>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                              ((expr_< table_<unspecified_<X>,SX>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                              ((expr_< table_<unspecified_<B>,SB>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                              ((expr_< table_<unspecified_<C>,SC>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                              ((expr_< table_<unspecified_<D>,SD>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                              )
  {
    typedef nt2::solve_lse_ip_return<A, B, C, D> result_type; 
    BOOST_DISPATCH_FORCE_INLINE result_type operator()(A& a, X&x, B const & b, C& c, D & d) const
    {
      return nt2::solve_lse_ip_return<A,B, C, D>(a,b,c,d);
    }
  };  
} }

namespace nt2 
{
  //============================================================================
  // lse actual functor : precompute
  //============================================================================
  template<class A, class B, class C, class D> struct solve_lse_ip_return
  {
    typedef typename A::value_type                   type_t;
    typedef typename A::index_type                  index_t; 
    typedef typename meta::as_real<type_t>::type    rtype_t; 
    typedef nt2::table<type_t,nt2::matlab_index_>    ftab_t;
    typedef nt2::table<rtype_t,nt2::matlab_index_>  fbtab_t;
    typedef nt2::table<la_int,nt2::matlab_index_>   fitab_t;
    typedef nt2::table<type_t,index_t>                tab_t;
    typedef nt2::table<rtype_t,index_t>              rtab_t;
    typedef nt2::table<la_int,index_t>               itab_t;

    ////////////////////////////////////////////////////////////////////////////
    // General LSE solver
    //  A is            M x N            may be rank-deficient
    //  X is or will be N x nrhs
    //  B is            M x nrhs    
    ////////////////////////////////////////////////////////////////////////////
    solve_lse_ip_return(A& a, X& x, const B& b)
    {
      la_int m = nt2::height(a);
      la_int n = nt2::width(a);
      la_int p = nt2::height(b);
      la_int lda = nt2::leading_size(a);
      la_int ldb = nt2::leading_size(b);
      la_int info; 
      BOOST_ASSERT_MSG( (n == nt2::width(b)),
                          "In lse calls the number of columns of a must match the number of columns of b");
        BOOST_ASSERT_MSG( (n >= p),
                          "In lse calls the number of columns of a must be greater or equal to the number of rows of b");
        BOOST_ASSERT_MSG( (n == nt2::numel(x)),
                          "In lse calls the number of columns of a must match the number of elements of x)");
        BOOST_ASSERT_MSG( (n <= m+p),
                          "In lse calls n <= m+p");
        BOOST_ASSERT_MSG( (p == nt2::numel(d)),
                          "In lse calls the number of rows of b must match the number of elements of d");
        //      std::cout << "conventional" << std::endl;
        nt2::details::gglse(&m, &n, &p, a.begin(), &lda, b.begin(), &ldb, c.begin(), d.begin(), x.begin(), &info);
    }
    ~solve_lse_ip_return(){}
    la_int get_info()    const { return info; }
    tab_t  get_x()       const { return x;    }
  private:
    la_int            info;
    tab_t                x; 
    
  };
} 

#endif
