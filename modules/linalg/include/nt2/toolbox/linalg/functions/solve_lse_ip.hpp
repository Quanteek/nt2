/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2009-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_SOLVE_LSE_IP_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_SOLVE_LSE_IP_HPP_INCLUDED
#include <nt2/include/functor.hpp>

/*!
 * \ingroup algebra
 * \defgroup algebra_solve_lse_ip solve_lse_ip
 *
 * \par Description
 * Elementary Least square
 *
 * \par Header file
 * 
 * \code
 * #include <nt2/include/functions/solve_lse_ip.hpp>
 * \endcode
 * 
 * 
 * \synopsis
 *
 * \param a the matrix a on entry, destroyed on exit
 *
 * \param x solution on exit
 * 
 * \param b the second member(s) on entry (remain untouched)
 *
 * \par Notes
 *   Call the dedicated lapack routines available on the target.
 * \par
 *  
**/
//==============================================================================
// solve_lse_ip actual class forward declaration
//==============================================================================
namespace nt2 
{
  template<class A, class B, class C, class D> struct solve_lse_ip_return;
} 

namespace nt2 { namespace tag
  {         
    /*!
     * \brief Define the tag solve_lse_ip_ of functor solve_lse_ip
     *        in namespace nt2::tag for toolbox algebra
    **/
    struct solve_lse_ip_ : ext::unspecified_<solve_lse_ip_> { typedef ext::unspecified_<solve_lse_ip_> parent; };
  }

  template<class A, class B, class C, class D>
  BOOST_FORCEINLINE nt2::solve_lse_ip_return<A, X, B> solve_lse_ip(A &a, B &b,  C &c, D &d)
  {
    return typename nt2::make_functor<tag::solve_lse_ip_, A>::type()(a,b,c,d);
  }

}

#endif

