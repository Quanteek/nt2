/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_WILKINSON_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_WILKINSON_HPP_INCLUDED
#include <nt2/include/functor.hpp>

/*! 
 * \ingroup algebra
 * \defgroup algebra_wilkinson wilkinson
 *
 * \par Description
 * compute a wilkinson matrix
 * wilkinson(n) is J. H. Wilkinson's eigenvalue test matrix, Wn+.
 * It is a symmetric, tridiagonal matrix with pairs of nearly,
 * but not exactly, equal eigenvalues.
 * The most frequently used case is wilkinson(21).
 * for example, wilkinson(7) is
 * 
 * 3  1  0  0  0  0  0
 * 1  2  1  0  0  0  0
 * 0  1  1  1  0  0  0
 * 0  0  1  0  1  0  0
 * 0  0  0  1  1  1  0
 * 0  0  0  0  1  2  1
 * 0  0  0  0  0  1  3
 *
 * \par Header file
 * 
 * \code
 * #include <nt2/include/functions/wilkinson.hpp>
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
// wilkinson actual class forward declaration
//==============================================================================

namespace nt2 { namespace tag
  {         
    /*!
     * \brief Define the tag wilkinson_ of functor wilkinson
     *        in namespace nt2::tag for toolbox algebra
    **/
    struct wilkinson_ : ext::unspecified_<wilkinson_> { typedef ext::unspecified_<wilkinson_> parent; };
  }
  
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION(tag::wilkinson_, wilkinson, 1)  

}

#endif

