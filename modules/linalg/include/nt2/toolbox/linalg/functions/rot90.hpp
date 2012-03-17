/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2009-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_ROT90_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_ROT90_HPP_INCLUDED
#include <nt2/include/functor.hpp>

/*!
 * \ingroup algebra
 * \defgroup algebra_rot90 rot90
 *
 * \par Description
 * Elementary Least square
 *
 * \par Header file
 * 
 * \code
 * #include <nt2/include/functions/rot90.hpp>
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
// rot90 actual class forward declaration
//==============================================================================

namespace nt2 { namespace tag
  {         
    /*!
     * \brief Define the tag rot90_ of functor rot90
     *        in namespace nt2::tag for toolbox algebra
    **/
    struct rot90_ : ext::unspecified_<rot90_> { typedef ext::unspecified_<rot90_> parent; };
  }
  
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION(tag::rot90_, rot90, 2)

}

#endif

