/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_CLEMENT_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_CLEMENT_HPP_INCLUDED
#include <nt2/include/functor.hpp>

/*! 
 * \ingroup algebra
 * \defgroup algebra_clement clement
 *
 * \par description
 *   clement(n,k) is a tridiagonal matrix with zero diagonal 
 *   entries and known eigenvalues. it is singular if n is odd. about 64
 *   percent of the entries of the inverse are zero. the eigenvalues
 *   are plus and minus the numbers n-1, n-3, n-5, ..., (1 or 0).
 *   for k = 0 (the default) the matrix is unsymmetric, while for
 *   k = 1 it is symmetric. clement(n,1) is diagonally similar 
 *   to clement(n).
 *
 *   note:
 *   similar properties hold for tridiag(x,y,z) where 
 *   y = zeros(n,1). the eigenvalues still come in plus/minus pairs but 
 *   they are not known explicitly.
 *
 *   references:
 *   [1] p.a. clement, a class of triple-diagonal matrices for test
 *   purposes, siam review, 1 (1959), pp. 50-52.
 *   [2] o. taussky and j. todd, another look at a matrix of mark kac,
 *   linear algebra and appl., 150 (1991), pp. 341-360.
 *
 *
 * \par Header file
 * 
 * \code
 * #include <nt2/include/functions/clement.hpp>
 * \endcode
 * 
 * 
 * \synopsis
 *
 * \param n order of the matrix output
 * \param k =  0 or 1
 * 
 *  
**/
//==============================================================================
// clement actual class forward declaration
//==============================================================================

namespace nt2 { namespace tag
  {         
    /*!
     * \brief Define the tag clement_ of functor clement
     *        in namespace nt2::tag for toolbox algebra
    **/
    struct clement_ : ext::unspecified_<clement_> { typedef ext::unspecified_<clement_> parent; };
  }
  
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION(tag::clement_, clement, 2)

}

#endif

