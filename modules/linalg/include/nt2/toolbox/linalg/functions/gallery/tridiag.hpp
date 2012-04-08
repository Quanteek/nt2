/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_TRIDIAG_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_TRIDIAG_HPP_INCLUDED
#include <nt2/include/functor.hpp>

/*! 
 * \ingroup algebra
 * \defgroup algebra_tridiag tridiag
 *
 * \par Description
 * compute a tridiag matrix expression from the three diagonals
 *
 * \par Header file
 * 
 * \code
 * #include <nt2/include/functions/tridiag.hpp>
 * \endcode
 * 
 * 
 * \synopsis
 *
 * \param c, d, e superdiag,  diag,  underdiag
 * 
 *  
**/
//==============================================================================
// tridiag actual class forward declaration
//==============================================================================

namespace nt2 { namespace tag
  {         
    /*!
     * \brief Define the tag tridiag_ of functor tridiag
     *        in namespace nt2::tag for toolbox algebra
    **/
    struct tridiag_ : ext::unspecified_<tridiag_> { typedef ext::unspecified_<tridiag_> parent; };
  }
  
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION(tag::tridiag_, tridiag, 3)  

}

#endif

