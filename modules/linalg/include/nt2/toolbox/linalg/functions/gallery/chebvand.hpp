/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_CHEBVAND_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_CHEBVAND_HPP_INCLUDED
#include <nt2/include/functor.hpp>

/*! 
 * \ingroup algebra
 * \defgroup algebra_chebvand chebvand
 *
 * \par Description
 * chebvand vandermonde-like matrix for the chebyshev polynomials.
 *    c = chebvand(p), where p is a vector, produces the
 *    (primal) chebyshev vandermonde matrix based on the points p:
 *       c(i,j) = t_{i-1}(p(j)), where t_{i-1} is the chebyshev
 *       polynomial of degree i-1.
 *    chebvand(m,p) is a rectangular version of
 *    chebvand(p) with m rows.
 *    special case: if p is a scalar, then p equally spaced points on
 *       [0,1] are used.
      
 *    reference:
 *    [1] n. j. higham, stability analysis of algorithms for solving confluent
 *        vandermonde-like systems, siam j. matrix anal. appl., 11 (1990),
 *        pp. 23-41.
 *
 * \par Header file
 * 
 * \code
 * #include <nt2/include/functions/chebvand.hpp>
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
// chebvand actual class forward declaration
//==============================================================================

namespace nt2 { namespace tag
  {         
    /*!
     * \brief Define the tag chebvand_ of functor chebvand
     *        in namespace nt2::tag for toolbox algebra
    **/
    struct chebvand_ : ext::unspecified_<chebvand_> { typedef ext::unspecified_<chebvand_> parent; };
  }
  
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION(tag::chebvand_, chebvand, 1)
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION(tag::chebvand_, chebvand, 2)

}

#endif

