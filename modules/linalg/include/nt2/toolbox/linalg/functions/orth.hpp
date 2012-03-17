/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2009-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_ORTH_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_ORTH_HPP_INCLUDED
#include <nt2/include/functor.hpp>

/*!
 * \ingroup algebra
 * \defgroup algebra_orth orth
 *
 * \par Description
 * Elementary Least square
 *
 * \par Header file
 * 
 * \code
 * #include <nt2/include/functions/orth.hpp>
 * \endcode
 * 
 * 
 * \synopsis
 *
 * \param a the matrix a on entry, destroyed on exit
 *
 * \param b the second member(s) b on entry, solution on exit
 * 
 * \par Notes
 *   Call the dedicated lapack routines available on the target.
 * \par
 *  
**/
//==============================================================================
// orth actual class forward declaration
//==============================================================================

namespace nt2 { namespace tag
  {         
    /*!
     * \brief Define the tag orth_ of functor orth
     *        in namespace nt2::tag for toolbox algebra
    **/
    struct orth_ : ext::unspecified_<orth_> { typedef ext::unspecified_<orth_> parent; };
  }
  
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION(tag::orth_, orth, 1)
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION(tag::orth_, orth, 2)

}

#endif

