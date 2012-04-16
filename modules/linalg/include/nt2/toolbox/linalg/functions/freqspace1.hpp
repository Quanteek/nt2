/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2009-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_FREQSPACE1_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_FREQSPACE1_HPP_INCLUDED
#include <nt2/include/functor.hpp>

/*!
 * \ingroup algebra
 * \defgroup algebra_freqspace1 freqspace1
 *
 * \par Description
 * Frequency spacing for frequency response. 1D case
 *
 * \par Header file
 * 
 * \code
 * #include <nt2/include/functions/freqspace1.hpp>
 * \endcode
 * 
 *   f =  freqspace1(n, _, as<T>())
 *   f =  freqspace1(n, as<T>())
 *
 *   f = freqspace1(n, as<T>()) returns the 1-d frequency vector f assuming n
 *                              equally spaced points around the unit circle.  
 *                              f = _(0, 2/n, 1).
 *   f = freqspace1(n, _, as<T>()) returns all n equally spaced points.
 *                                 In this case, f = _(0, 2/n, 2*(n-1)/n).
 *
 * T can be any floating type
**/
//==============================================================================
// freqspace1 actual class forward declaration
//==============================================================================

namespace nt2 { namespace tag
  {         
    /*!
     * \brief Define the tag freqspace1_ of functor freqspace1
     *        in namespace nt2::tag for toolbox algebra
    **/
    struct freqspace1_ : ext::unspecified_<freqspace1_> { typedef ext::unspecified_<freqspace1_> parent; };
  }
  
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION(tag::freqspace1_, freqspace1, 2)
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION(tag::freqspace1_, freqspace1, 3)

}

#endif

