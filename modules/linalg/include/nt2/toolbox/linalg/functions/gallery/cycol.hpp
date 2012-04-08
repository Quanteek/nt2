/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_CYCOL_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_CYCOL_HPP_INCLUDED
#include <nt2/include/functor.hpp>

/*! 
 * \ingroup algebra
 * \defgroup algebra_cycol cycol
 *
 * \par Description
 * compute a cycol matrix
 * cycol  matrix whose columns repeat cyclically.
 *    a = cycol(m, n, k) is an m-by-n matrix of the form 
 *    a = b(1:m,1:n) where b = [c c c...] and c = randn(m, k). thus a's 
 *    columns repeat cyclically, and a has rank at most k. k need not 
 *    divide n. k defaults to round(n/4).
 *    cycol(n,k), where n is a scalar, is the same as 
 *    cycol(n, n, k).
      
 *    note:
 *    this type of matrix can lead to underflow problems for gaussian
 *    elimination. see the reference below.
 * 
 *    reference:
 *    [1] na digest volume 89, issue 3 (january 22, 1989).
 * 
 * Parameter n specifies dimension: m-by-n.
 * \par Header file
 * 
 * \code
 * #include <nt2/include/functions/cycol.hpp>
 * \endcode
 * 
 * 
 * \synopsis
 *
 * \param n, m order of the matrix output
 * 
 *  
**/
//==============================================================================
// cycol actual class forward declaration
//==============================================================================

namespace nt2 { namespace tag
  {         
    /*!
     * \brief Define the tag cycol_ of functor cycol
     *        in namespace nt2::tag for toolbox algebra
    **/
    struct cycol_ : ext::unspecified_<cycol_> { typedef ext::unspecified_<cycol_> parent; };
  }
  

  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION(tag::cycol_, cycol, 4)
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION(tag::cycol_, cycol, 3)  

}

#endif

