/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_DORR_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_DORR_HPP_INCLUDED
#include <nt2/include/functor.hpp>

/*! 
 * \ingroup algebra
 * \defgroup algebra_dorr dorr
 *
 * \par description
 * dorr dorr matrix (sparse).
 *    dorr(n,theta, c, d, e) returns the vectors defining a row
 *    diagonally dominant, tridiagonal n-by-n matrix that is ill-conditioned
 *    for small values of theta >= 0. theta defaults to 0.01.
 * 
 *    a = dorr(n,theta) returns the dorr matrix itself.  
 *    this is the same as
 *        dorr(n,theta, c, d, e);
 *        a = tridiag(c,d,e);
 * 
 *    the columns of inv(a) vary greatly in norm.
 *    the amount of diagonal dominance, ignoring rounding errors, is:
 *         compar(c)*ones(n,1) = theta*sqr(n+1)*catv(1, zeros(1, n-2), 1).
      
 *    reference:
 *    [1] f. w. dorr, an example of ill-conditioning in the numerical
 *        solution of singular perturbation problems, math. comp.,
 *        25 (1971), pp. 271-283.
 *
 * \par header file
 * 
 * \code
 * #include <nt2/include/functions/dorr.hpp>
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
// dorr actual class forward declaration
//==============================================================================

namespace nt2 { namespace tag
  {         
    /*!
     * \brief Define the tag dorr_ of functor dorr
     *        in namespace nt2::tag for toolbox algebra
    **/
    struct dorr_ : ext::unspecified_<dorr_> { typedef ext::unspecified_<dorr_> parent; };
  }
  
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION(tag::dorr_, dorr, 5)  
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION(tag::dorr_, dorr, 2)  

}

#endif

