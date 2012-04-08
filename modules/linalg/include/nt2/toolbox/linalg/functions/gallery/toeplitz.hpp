/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_TOEPLITZ_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_TOEPLITZ_HPP_INCLUDED
#include <nt2/include/functor.hpp>

/*! 
 * \ingroup algebra
 * \defgroup algebra_toeplitz toeplitz
 *
 * \par Description
 *   toeplitz(c,r) is a non-symmetric toeplitz matrix having c as its
 *   first column and r as its first row.
 *   the first c element wins against the first r element
 *
 *   toeplitz(r) is a symmetric (or hermitian) toeplitz matrix.
 *
 *   see also hankel.
 *
 * \par Header file
 * 
 * \code
 * #include <nt2/include/functions/toeplitz.hpp>
 * \endcode
 * 
 * 
 * \synopsis
 *
 * \param n order of the matrix output
 * 
 *  
**/
//==============================================================================
// toeplitz actual class forward declaration
//==============================================================================

namespace nt2 { namespace tag
  {         
    /*!
     * \brief Define the tag toeplitz_ of functor toeplitz
     *        in namespace nt2::tag for toolbox algebra
    **/
    struct toeplitz_ : ext::unspecified_<toeplitz_> { typedef ext::unspecified_<toeplitz_> parent; };
  }
  
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION(tag::toeplitz_, toeplitz, 1)
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION(tag::toeplitz_, toeplitz, 2)

}

#endif

