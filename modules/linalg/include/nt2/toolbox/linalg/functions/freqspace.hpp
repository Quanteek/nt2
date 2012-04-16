/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2009-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_FREQSPACE_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_FREQSPACE_HPP_INCLUDED
#include <nt2/include/functor.hpp>

/*!
 * \ingroup algebra
 * \defgroup algebra_freqspace freqspace
 *
 * \par Description
 * Elementary Least square
 *
 * \par Header file
 * 
 * \code
 * #include <nt2/include/functions/freqspace.hpp>
 * \endcode
 * 
 * 
 * \synopsis
 * in matlab freqspace return a result that depend of the output type
 * I have written two different functors freqspace1 and freqspace2
 * TODO combine them 
 *  
**/
//==============================================================================
// freqspace actual class forward declaration
//==============================================================================

namespace nt2 { namespace tag
  {         
    /*!
     * \brief Define the tag freqspace_ of functor freqspace
     *        in namespace nt2::tag for toolbox algebra
    **/
    struct freqspace_ : ext::unspecified_<freqspace_> { typedef ext::unspecified_<freqspace_> parent; };
  }
  
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION(tag::freqspace_, freqspace, 2)
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION(tag::freqspace_, freqspace, 3)
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION(tag::freqspace_, freqspace, 4)

}

#endif

