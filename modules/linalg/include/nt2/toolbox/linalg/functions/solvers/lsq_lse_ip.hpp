/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_SOLVERS_LSQ_LSE_IP_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_SOLVERS_LSQ_LSE_IP_HPP_INCLUDED

#include <nt2/include/functions/lsq_lse_ip.hpp>
#include <nt2/include/functions/leading_size.hpp>
#include <nt2/toolbox/linalg/details/lapack/gglse.hpp>
#include <nt2/table.hpp>
#include <nt2/include/functions/max.hpp>
#include <nt2/include/functions/min.hpp>
#include <nt2/include/functions/areofsameheight.hpp>
#include <nt2/include/functions/height.hpp>
#include <nt2/include/functions/width.hpp>


//==============================================================================
// lse actual functor forward declaration
//==============================================================================
namespace nt2
{
  template<class A, class B, class C, class D, class X> struct lsq_lse_ip_return;
} 

namespace nt2 { namespace ext
{
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::lsq_lse_ip_, tag::cpu_, 
                              (A)(SA)(X)(SX)(B)(SB)(C)(SC)(D)(SD(X(SX)),
                              ((expr_< table_<unspecified_<A>,SA>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                              ((expr_< table_<unspecified_<X>,SX>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                              ((expr_< table_<unspecified_<B>,SB>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                              ((expr_< table_<unspecified_<C>,SC>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                              ((expr_< table_<unspecified_<D>,SD>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                              ((expr_< table_<unspecified_<X>,SX>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                             )
  {
    typedef nt2::lsq_lse_ip_return<A, B, C, D, X> result_type; 
    BOOST_DISPATCH_FORCE_INLINE result_type operator()(A& a, B const & b, C& c, D & d, X& x) const
    {
      return nt2::lsq_lse_ip_return<A,B, C, D, X>(a,b,c,d,x);
    }
  };  
} }

namespace nt2 
{
  //============================================================================
  // lse actual functor : precompute
  //============================================================================
  template<class A, class B, class C, class D, class X> struct lsq_lse_ip_return
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
    lsq_lse_ip_return(A& a, B& b, C&c, D&d, X&x)
    {
      la_int m = nt2::height(a);
      la_int n = nt2::width(a);
      la_int p = nt2::height(b);
      la_int lda = nt2::leading_size(a);
      la_int ldb = nt2::leading_size(b);
      BOOST_ASSERT_MSG( (n == nt2::width(b)),"In lse calls the number of columns of a must match the number of columns of b");
      BOOST_ASSERT_MSG( (n >= p),"In lse calls the number of columns of a must be greater or equal to the number of rows of b");
      BOOST_ASSERT_MSG( (n == nt2::numel(x)),"In lse calls the number of columns of a must match the number of elements of x)");
      BOOST_ASSERT_MSG( (n <= m+p),"In lse calls n <= m+p");
      BOOST_ASSERT_MSG( (p == nt2::numel(d)),"In lse calls the number of rows of b must match the number of elements of d");
      nt2::details::gglse(&m, &n, &p, a.raw(), &lda, b.raw(), &ldb, c.raw(), d.raw(), x.raw(), &info);
    }
    ~lsq_lse_ip_return(){}
    la_int get_info()    const { return info; }
  private:
    la_int            info;
  };
} 

#endif
