/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_CIRCUL_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_CIRCUL_HPP_INCLUDED
#include <nt2/include/functor.hpp>

/*! 
 * \ingroup algebra
 * \defgroup algebra_circul circul
 *
 * \par Description
 * compute a circul matrix
 *   c = circul(v) is the circulant matrix whose first row is v.
 *   a circulant matrix has the property that each row is obtained
 *   from the previous one by cyclically permuting the entries one step
 *   forward. it is a special toeplitz matrix in which the diagonals
 *   "wrap round". if v is a scalar, then c = circul(1:v).
 *
 *   the eigensystem of c (n-by-n) is known explicitly. if t is an nth
 *   root of unity, then the inner product of v with w = [1 t t^2 ... t^(n-1)]
 *   is an eigenvalue of c, and w(n:-1:1) is an eigenvector of c.
      
 *   reference:
 *   [1] p.j. davis, circulant matrices, john wiley, 1977.
 *
 * \par Header file
 * 
 * \code
 * #include <nt2/include/functions/circul.hpp>
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
// circul actual class forward declaration
//==============================================================================

namespace nt2 { namespace tag
  {         
    /*!
     * \brief Define the tag circul_ of functor circul
     *        in namespace nt2::tag for toolbox algebra
    **/
    struct circul_ : ext::unspecified_<circul_> { typedef ext::unspecified_<circul_> parent; };
  }
  
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION(tag::circul_, circul, 1)

}

#endif

