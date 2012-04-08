/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_CHOW_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_CHOW_HPP_INCLUDED
#include <nt2/include/functor.hpp>

/*! 
 * \ingroup algebra
 * \defgroup algebra_chow chow
 *
 * \par description
 *  a = chow(n,alpha,delta) returns a such that
 *  a = h(alpha) + delta*eye, where h(i,j) = alpha^(i-j+1).
 *  h(alpha) has p = floor(n/2) zero eigenvalues, the rest being
 *  4*alpha*cos( k*pi/(n+2) )^2, k=1:n-p.
 *  defaults: alpha = 1, delta = 0.
   
 *  references:
 *  [1] t. s. chow, a class of hessenberg matrices with known eigenvalues
 *      and inverses, siam review, 11 (1969), pp. 391-395.
 *  [2] g. fairweather, on the eigenvalues and eigenvectors of a class of
 *      hessenberg matrices, siam review, 13 (1971), pp. 220-221.
 *  [3] i. singh, g. poole and t. boullion, a class of hessenberg matrices
 *      with known pseudoinverse and drazin inverse, math. comp., 29 (1975),
 *      pp. 615-619.
 * \par Header file
 * 
 * \code
 * #include <nt2/include/functions/chow.hpp>
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
// chow actual class forward declaration
//==============================================================================

namespace nt2 { namespace tag
  {         
    /*!
     * \brief Define the tag chow_ of functor chow
     *        in namespace nt2::tag for toolbox algebra
    **/
    struct chow_ : ext::unspecified_<chow_> { typedef ext::unspecified_<chow_> parent; };
  }
  
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION(tag::chow_, chow, 1)
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION(tag::chow_, chow, 2)
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION(tag::chow_, chow, 3) 

}

#endif

