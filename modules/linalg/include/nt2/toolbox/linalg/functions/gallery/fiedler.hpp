/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_FIEDLER_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_FIEDLER_HPP_INCLUDED
#include <nt2/include/functor.hpp>

/*! 
 * \ingroup algebra
 * \defgroup algebra_fiedler fiedler
 *
 * \par Description
 *    a = fiedler(c), where c is an n-vector, is the n-by-n 
 *    symmetric matrix with elements abs(c(i)-c(j)). if c is a scalar, 
 *    then a = fiedler(_(1, c)).
 * 
 *    a has a dominant positive eigenvalue and all the other eigenvalues
 *    are negative. (Szego 1936)
 * 
 *    note: explicit formulas for inv(a) and det(a) are given in (todd 1977)
 *    and attributed to fiedler. these indicate that inv(a) is 
 *    tridiagonal except for nonzero (1,n) and (n,1) elements.
 *     
 *    References:
 *      [1] G. Szego, Solution to problem 3705, Amer. Math. Monthly,
 *        43 (1936), pp. 246-259.
 *      [2] J. Todd, Basic Numerical Mathematics, Vol. 2: Numerical Algebra,
 *        Birkhauser, Basel, and Academic Press, New York, 1977, p. 159.
 * 
 * \par Header file
 * 
 * \code
 * #include <nt2/include/functions/fiedler.hpp>
 * \endcode
 * 
 * 
 * \synopsis
 *
 * \param c n-vector or scalar integer
 * 
 *  
**/
//==============================================================================
// fiedler actual class forward declaration
//==============================================================================

namespace nt2 { namespace tag
  {         
    /*!
     * \brief Define the tag fiedler_ of functor fiedler
     *        in namespace nt2::tag for toolbox algebra
    **/
    struct fiedler_ : ext::unspecified_<fiedler_> { typedef ext::unspecified_<fiedler_> parent; };
  }
  
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION(tag::fiedler_, fiedler, 1)  

}

#endif

