/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_SOLVERS_SOLVE_SVD_IP_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_SOLVERS_SOLVE_SVD_IP_HPP_INCLUDED

#include <nt2/include/functions/solve_svd_ip.hpp>
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
// svd actual functor forward declaration
//==============================================================================
namespace nt2
{
  template<class A, class X, class B> struct solve_svd_ip_return;
} 

namespace nt2 { namespace ext
{
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::solve_svd_ip_, tag::cpu_, 
                              (A)(SA)(X)(SX)(B)(SB),
                              ((expr_< table_<unspecified_<A>,SA>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                              ((expr_< table_<unspecified_<X>,SX>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                              ((expr_< table_<unspecified_<B>,SB>,nt2::tag::terminal_,boost::mpl::long_<0> >))                              )
  {
    typedef nt2::solve_svd_ip_return<A, X, B> result_type; 
    BOOST_DISPATCH_FORCE_INLINE result_type operator()(A& a, X&x, B const & b ) const
    {
      return nt2::solve_svd_ip_return<A, X, B>(a,x,b);
    }
  };  
} }

namespace nt2 
{
  //============================================================================
  // svd actual functor : precompute
  //============================================================================
  template<class A, class B> struct solve_svd_ip_return
  {
    typedef long int                                 la_int; 
    typedef typename A::value_type                   type_t;
    typedef typename A::index_type                  index_t; 
    typedef typename meta::as_real<type_t>::type    rtype_t; 
    typedef nt2::table<type_t,nt2::matlab_index_>    ftab_t;
    typedef nt2::table<btype_t,nt2::matlab_index_>  fbtab_t;
    typedef nt2::table<la_int,nt2::matlab_index_>   fitab_t;
    typedef nt2::table<type_t,index_t>                tab_t;
    typedef nt2::table<rtype_t,index_t>              rtab_t;
    typedef nt2::table<la_int,index_t>               itab_t;

    ////////////////////////////////////////////////////////////////////////////
    // General SVD solver
    //  A is            M x N            may be rank-deficient
    //  X is or will be N x nrhs
    //  B is            M x nrhs    
    ////////////////////////////////////////////////////////////////////////////
    solve_svd_ip_return(A& a, X& x, const B& b)
    {
      const la_int ml = size(a, 1);
      const la_int nl = size(a, 2);
      const la_int nrhs = size(b, 2);
      const la_int lda = leading_size(a); 
      const la_int ldb = leading_size(b); 
      rtab_t s(of_size(nt2::min(ml, nl), 1)); 
      
      // typically is a non-square, so we need to create tmp x because is 
      //  x is n x nrhs, while b is m x nrhs.  we need to make copies of
      //  these so that the routine won't corrupt data around x and b
      
      if (ml != nl)
        {
          la_int mm =  std::max(std::max(ml,nl),1l);
          tab_t xtmp = b; //nt2::expand(b, nt2::of_size(mm, nrhs));
          nt2::details::gelsd(&ml, &nl, &nrhs, a.raw(), &lda, xtmp.raw(), &ldb,
                              s.raw(), &rcond, &rank, &info);                   
          x = xtmp; //(range(1, nl), range(1, nrhs));
          boost_assert_msg(info!= 0, "lapack error : gelsd in solve_svd_ip(1)");
        }
      else
        {
          x = b; 
          nt2::details::gelsd(&ml, &nl, &nrhs, a.raw(), &lda, x.raw(), &ldb,
                              s.raw(), &rcond, &rank, &info);                   
          boost_assert_msg(info!= 0, "lapack error : gelsd in solve_svd_ip(2)");
        }
    }
    ~solve_svd_ip_return(){}
    rtype_t get_rcond()  const { return rcond;}   
    la_int get_rank()    const { return rank; }
    la_int get_info()    const { return info; }
  private:
    la_int                             info;
    la_int                             rank; 
    rtype_t                           rcond;                                      
    
  };
} 

#endif
