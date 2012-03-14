/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2009-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_SOLVE_SVD_IP_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_SOLVE_SVD_IP_HPP_INCLUDED
#include <nt2/include/functor.hpp>

/*!
 * \ingroup algebra
 * \defgroup algebra_solve_svd_ip solve_svd_ip
 *
 * \par Description
 * Elementary Least square
 *
 * \par Header file
 * 
 * \code
 * #include <nt2/include/functions/solve_svd_ip.hpp>
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
// solve_svd_ip actual class forward declaration
//==============================================================================
namespace nt2 
{
  template<class A, class X, class B> struct solve_svd_ip_return;
} 

namespace nt2 { namespace tag
  {         
    /*!
     * \brief Define the tag solve_svd_ip_ of functor solve_svd_ip
     *        in namespace nt2::tag for toolbox algebra
    **/
    struct solve_svd_ip_ : ext::unspecified_<solve_svd_ip_> { typedef ext::unspecified_<solve_svd_ip_> parent; };
  }

  template<class A, class X, class B>
  BOOST_FORCEINLINE nt2::solve_svd_ip_return<A, X, B> solve_svd_ip(A &a, X &x, B &b)
  {
    return typename nt2::make_functor<tag::solve_svd_ip_, A>::type()(a, x, b);
  }

}

#endif

