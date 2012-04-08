/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_CHEBSPEC_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_CHEBSPEC_HPP_INCLUDED
#include <nt2/include/functor.hpp>

/*! 
 * \ingroup algebra
 * \defgroup algebra_chebspec chebspec
 *
 * \par Description
 * chebspec chebyshev spectral differentiation matrix.
 *    c = chebspec(n,k) is a chebyshev spectral
 *       differentiation matrix of order n.  k = 0 (the default) or 1.
 *    for k = 0 ("no boundary conditions"), c is nilpotent, with
 *       c^n = 0 and it has the null vector ones(n,1).
 *       c is similar to a jordan block of size n with eigenvalue zero.
 *    for k = 1, c is nonsingular and well conditioned, and its
 *       eigenvalues have negative real parts.
 *    for both k, the computed eigenvector matrix x from eig is
 *       ill-conditioned (mesh(real(x)) is interesting).
      
 *    references:
 *    [1] c. canuto, m. y. hussaini, a. quarteroni and t. a. zang,
 *        spectral methods in fluid dynamics, springer-verlag, berlin,
 *        1988, p. 69.
 *    [2] l. n. trefethen and m. r. trummer, an instability phenomenon in
 *        spectral methods, siam j. numer. anal., 24 (1987), pp. 1008-1023.
 *    [3] d. funaro, computing the inverse of the chebyshev collocation
 *        derivative, siam j. sci. stat. comput., 9 (1988), pp. 1050-1057.
 *
 * \par Header file
 * 
 * \code
 * #include <nt2/include/functions/chebspec.hpp>
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
// chebspec actual class forward declaration
//==============================================================================

namespace nt2 { namespace tag
  {         
    /*!
     * \brief Define the tag chebspec_ of functor chebspec
     *        in namespace nt2::tag for toolbox algebra
    **/
    struct chebspec_ : ext::unspecified_<chebspec_> { typedef ext::unspecified_<chebspec_> parent; };
  }
  
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION(tag::chebspec_, chebspec, 1)
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION(tag::chebspec_, chebspec, 2)

}

#endif
