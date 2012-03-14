/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_SOLVERS_SOLVE_LU_IP_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_SOLVERS_SOLVE_LU_IP_HPP_INCLUDED

#include <nt2/include/functions/solve_lu_ip.hpp>
#include <nt2/include/functions/numel.hpp>
#include <nt2/include/functions/leading_size.hpp>
#include <nt2/toolbox/linalg/details/lapack/gesv.hpp>
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
#include <iostream>
#include <vector>


//==============================================================================
// svd actual functor forward declaration
//==============================================================================
namespace nt2
{
  template<class A, class B> struct solve_lu_ip_return;
} 

namespace nt2 { namespace ext
{
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::solve_lu_ip_, tag::cpu_, 
                              (A)(SA)(B)(SB),
                              ((expr_< table_<unspecified_<A>,SA>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                              ((expr_< table_<unspecified_<B>,SB>,nt2::tag::terminal_,boost::mpl::long_<0> >))   
                              )
  {
    typedef nt2::solve_lu_ip_return<A, B> result_type; 
    BOOST_DISPATCH_FORCE_INLINE result_type operator()(A& a,B& b ) const
    {
      return nt2::solve_lu_ip_return<A, B>(a, b);
    }
  };  
} }

namespace nt2 
{
  //============================================================================
  // svd actual functor : precompute
  //============================================================================
  template<class A, class B> struct solve_lu_ip_return
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
    // General LU Solver
    //  A is            N x N
    //  B is            N x nrhs
    ////////////////////////////////////////////////////////////////////////////
    solve_lu_ip_return(A& a, B& b) : ipiv(nt2::of_size(height(a), 1))
    {
      //       std::cout << nt2::type_id<itab_t>() << std::endl; 
      BOOST_ASSERT_MSG(nt2::issquare(a), "matrix A is not square");
      BOOST_ASSERT_MSG(nt2::areofsameheight(a, b), "A and X have different heights");
      la_int info;
      la_int Ml   = nt2::height(a);
      la_int K    = nt2::width(b);
      la_int lda  = nt2::leading_size(a);
      la_int ldx  = nt2::leading_size(b);
      nt2::details::gesv (&Ml, &K, a.raw(), &lda, ipiv.raw(), b.raw(), &ldx, &info);
      BOOST_ASSERT_MSG(info == 0, "Lapack error : gesv in solve_lu_ip");
    }
    ~solve_lu_ip_return(){}
    itab_t get_piv()     const { return ipiv; }
    la_int get_info()    const { return info; }
  private:
    fitab_t                            ipiv;
    la_int                             info; 
    
  };
} 

#endif
