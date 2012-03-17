/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2009-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_SVD_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_SVD_HPP_INCLUDED

#include <nt2/include/functor.hpp>

/*!
 * \ingroup algebra
 * \defgroup algebra_svd svd
 *
 * \par Description
 * Elementary Least square
 *
 * \par Header file
 * 
 * \code
 * #include <nt2/include/functions/svd.hpp>
 * \endcode
 * 
 * 
 * \synopsis
 *
 * \param a the matrix a on entry, destroyed on exit
 *
 * \par Notes
 *   Call the dedicated lapack routines available on the target.
 * \par
 *  
**/
//==============================================================================
// svd actual class forward declaration
//==============================================================================
namespace nt2 
{
  template<class A> struct svd_return;
} 

namespace nt2 { namespace tag
  {         
    /*!
     * \brief Define the tag svd_ of functor svd
     *        in namespace nt2::tag for toolbox algebra
    **/
    struct svd_ : ext::unspecified_<svd_> { typedef ext::unspecified_<svd_> parent; };
  }

  template<class A>
  BOOST_FORCEINLINE nt2::svd_return<A> svd(A &a, const char & jobz = 'A')
  {
    return typename nt2::make_functor<tag::svd_, A>::type()(a, jobz);
  }

}

#endif

