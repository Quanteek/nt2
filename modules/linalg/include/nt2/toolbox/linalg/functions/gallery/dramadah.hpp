/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_DRAMADAH_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_DRAMADAH_HPP_INCLUDED
#include <nt2/include/functor.hpp>

/*! 
 * \ingroup algebra
 * \defgroup algebra_dramadah dramadah
 *
 * \par Description
 * dramadah matrix of zeros and ones with large determinant or inverse.
 *    a = dramadah(n,k) is an n-by-n (0,1) matrix for which
 *    mu(a) = norm(inv(a),'fro') or det(a) is relatively large.
 * 
 *    k = 1: (default)
 *       a is toeplitz, with abs(det(a)) = 1, and mu(a) > c(1.75)^n,
 *       where c is a constant. inv(a) has integer entries.
 *    k = 2:
 *       a is upper triangular and toeplitz. inv(a) has integer entries.
 *    k = 3:
 *       a has maximal determinant among (0,1) lower hessenberg matrices.
 *       det(a) = the n'th fibonacci number. a is toeplitz.
 *       the eigenvalues have an interesting distribution in the complex
 *       plane.
 * 
 *    an anti-hadamard matrix a is a matrix with elements 0 or 1 for
 *    which mu(a) = norm(inv(a),'fro') is maximal.  for k = 1,2 this function
 *    returns matrices with mu(a) relatively large, though not necessarily
 *    maximal.
      
 *    references:
 *    [1] r. l. graham and n. j. a. sloane, anti-hadamard matrices,
 *        linear algebra and appl., 62 (1984), pp. 113-137.
 *    [2] l. ching, the maximum determinant of an nxn lower hessenberg
 *        (0,1) matrix, linear algebra and appl., 183 (1993), pp. 147-153.
 * 
 *    nicholas j. higham, dec 1999.
 *
 * \par Header file
 * 
 * \code
 * #include <nt2/include/functions/dramadah.hpp>
 * \endcode
 * 
 * 
 * \synopsis
 *
 * \param n matrix order
 * \param k,  k =  1 or 2 or 3
 * 
 *  
**/
//==============================================================================
// dramadah actual class forward declaration
//==============================================================================

namespace nt2 { namespace tag
  {         
    /*!
     * \brief Define the tag dramadah_ of functor dramadah
     *        in namespace nt2::tag for toolbox algebra
    **/
    struct dramadah_ : ext::unspecified_<dramadah_> { typedef ext::unspecified_<dramadah_> parent; };
  }
  
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION(tag::dramadah_, dramadah, 1)  
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION(tag::dramadah_, dramadah, 2)  

}

#endif

