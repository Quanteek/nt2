/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_MAGIC_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_MAGIC_HPP_INCLUDED
#include <nt2/include/functor.hpp>

/*! 
 * \ingroup algebra
 * \defgroup algebra_magic magic
 *
 * \par Description
 *   magic  magic square.
 *   magic(n) is an n-by-n matrix constructed from the integers
 *   1 through n^2 with equal row, column, and diagonal sums.
 *   produces valid magic squares for all n > 0 except n = 2.
 *
 * \par Header file
 * 
 * \code
 * #include <nt2/include/functions/magic.hpp>
 * \endcode
 * 
 * 
 * \synopsis
 *
 * \param n order of the matrix output (n != 2)
 * 
 *  
**/
//==============================================================================
// magic actual class forward declaration
//==============================================================================

namespace nt2 { namespace tag
  {         
    /*!
     * \brief Define the tag magic_ of functor magic
     *        in namespace nt2::tag for toolbox algebra
    **/
    struct magic_ : ext::unspecified_<magic_> { typedef ext::unspecified_<magic_> parent; };
  }
  
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION(tag::magic_, magic, 1)

}

#endif

