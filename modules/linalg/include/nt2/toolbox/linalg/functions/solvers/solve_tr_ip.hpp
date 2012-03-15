/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_SOLVERS_SOLVE_TR_IP_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_SOLVERS_SOLVE_TR_IP_HPP_INCLUDED

#include <nt2/include/functions/solve_tr_ip.hpp>
#include <nt2/include/functions/numel.hpp>
#include <nt2/include/functions/leading_size.hpp>
#include <nt2/toolbox/linalg/details/lapack/trtrs.hpp>
#include <nt2/table.hpp>
#include <nt2/include/constants/zero.hpp>
#include <nt2/include/constants/mone.hpp>
#include <nt2/include/functions/eps.hpp>
#include <nt2/include/functions/rec.hpp>
#include <nt2/include/functions/max.hpp>
#include <nt2/include/functions/min.hpp>
#include <nt2/include/functions/issquare.hpp>
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


//==============================================================================
// svd actual functor forward declaration
//==============================================================================
namespace nt2
{
  template<class A, class B> struct solve_tr_ip_return;
} 

namespace nt2 { namespace ext
{
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::solve_tr_ip_, tag::cpu_, 
                              (A)(SA)(B)(SB)(C),
                              ((expr_< table_<unspecified_<A>,SA>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                              ((expr_< table_<unspecified_<B>,SB>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                              (scalar_<integer_<C> >)
                              (scalar_<integer_<C> >)
                              (scalar_<integer_<C> >)
                              )
  {
    typedef nt2::solve_tr_ip_return<A, B> result_type; 
    BOOST_DISPATCH_FORCE_INLINE result_type operator()(A& a,B& b,
                                                       const char & uplo,
                                                       const char & trans,
                                                       const char & diag ) const
    {
      return nt2::solve_tr_ip_return<A, B>(a, b);
    }
  };  
} }

namespace nt2 
{
  //============================================================================
  // svd actual functor : precompute
  //============================================================================
  template<class A, class B> struct solve_tr_ip_return
  {
    typedef typename A::value_type                   type_t;
    typedef typename A::index_type                  index_t; 
    typedef typename meta::as_real<type_t>::type    btype_t; 
    typedef nt2::table<type_t,nt2::matlab_index_>    ftab_t;
    typedef nt2::table<btype_t,nt2::matlab_index_>  fbtab_t;
    typedef nt2::table<la_int,nt2::matlab_index_>   fitab_t;
    typedef nt2::table<type_t,index_t>                tab_t;
    typedef nt2::table<btype_t,index_t>              btab_t;
    typedef nt2::table<la_int,index_t>               itab_t;
    
    ////////////////////////////////////////////////////////////////////////////
    //  Solve triangular systems
    //  A is            N x N
    //  X is            N x nrhs
    //  need sentimentally a triangular matrix. Only the lower
    //  triangular part is used by default (uplo = 'l' : use 'u' to get the
    //  upper part used)
    //  solve a*xx =  b, use trans = 't' to solve trans(a)*xx = b
    //  do not assume diag elements are special, use diag = 'u' to suppose ones
    //  on the diagonal and no access to a diagonal elements
    ////////////////////////////////////////////////////////////////////////////
    solve_tr_ip_return(A& a, B& bx,
                       const char & uplo = 'l',
                       const char& trans = 'n',      
                       const char& diag  = 'n' ) 
    {
      BOOST_ASSERT_MSG(nt2::issquare(a), "matrix a is not square");
      BOOST_ASSERT_MSG(nt2::areofsameheight(a, bx), "a and x have different heights");
      const la_int m = height(a);
      const la_int k = width(bx);
      const la_int lda = leading_size(a);
      const la_int ldx = leading_size(bx);
      nt2::details::trtrs(&uplo, &trans, &diag, &m, &k, a.raw(), &lda, bx.raw(), &ldx, &info);
      BOOST_ASSERT_MSG(info == 0, "lapack error : gels in solve_tr_ip");
    }
    ~solve_tr_ip_return(){}
    la_int get_info()    const { return info; }
  private:
    la_int  info; 
    
  };
} 

#endif
