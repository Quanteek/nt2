/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2009-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_RREF_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_RREF_HPP_INCLUDED
#include <nt2/include/functor.hpp>

/*!
 * \ingroup algebra
 * \defgroup algebra_rref rref
 *
 * \par Description
 * reduced echelon form
 *
 * \par Header file
 * 
 * \code
 * #include <nt2/include/functions/rref.hpp>
 * \endcode
 * 
 * Reduced row echelon form.
 *
 *   tie(r,jb) = rref(a)
 *   r is the reduced row echelon form of a.
 *   jb is a vector so that:
 *       r = length(jb) is this algorithm's idea of the rank of a,
 *       x(jb) are the bound variables in a linear system, ax = b,
 *       a(_,jb) is a basis for the range of a,
 *       r(_(1,r),jb) is the r-by-r identity matrix.
 *
 *   tie(r,jb) = rref(a,tol) uses the given tolerance in the rank tests.
 *
 *   roundoff errors may cause this algorithm to compute a different
 *   value for the rank than rank, orth and null.
/ * 
 * \synopsis
 * tie(r, jb) = rref(a, tol = Eps)
 *  
**/
//==============================================================================
// rref actual class forward declaration
//==============================================================================

namespace nt2 { namespace tag
  {         
    /*!
     * \brief Define the tag rref_ of functor rref
     *        in namespace nt2::tag for toolbox algebra
    **/
    struct rref_ : ext::unspecified_<rref_> { typedef ext::unspecified_<rref_> parent; };
  }
  
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION(tag::rref_, rref, 2)
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION(tag::rref_, rref, 1)

}

#endif

