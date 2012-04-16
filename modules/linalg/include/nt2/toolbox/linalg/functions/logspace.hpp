/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2009-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_LOGSPACE_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_LOGSPACE_HPP_INCLUDED
#include <nt2/include/functor.hpp>

/*!
 * \ingroup algebra
 * \defgroup algebra_logspace logspace
 *
 * \par Description
 * Elementary Least square
 *
 * \par Header file
 * 
 * \code
 * #include <nt2/include/functions/logspace.hpp>
 * \endcode
 * 
 * 
 * \synopsis
 * in matlab logspace return a result that depend of the output type
 * I have written two different functors logspace1 and logspace2
 * TODO combine them 
 *  
**/
//==============================================================================
// logspace actual class forward declaration
//==============================================================================

namespace nt2 { namespace tag
  {         
    /*!
     * \brief Define the tag logspace_ of functor logspace
     *        in namespace nt2::tag for toolbox algebra
    **/
    struct logspace_ : ext::unspecified_<logspace_> { typedef ext::unspecified_<logspace_> parent; };
  }
  
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION(tag::logspace_, logspace, 2)
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION(tag::logspace_, logspace, 3)

}

#endif

