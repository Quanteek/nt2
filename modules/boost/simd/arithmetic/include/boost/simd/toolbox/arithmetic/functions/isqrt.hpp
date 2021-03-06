//////////////////////////////////////////////////////////////////////////////
///   Copyright 2003 and onward LASMEA UMR 6602 CNRS/U.B.P Clermont-Ferrand
///   Copyright 2009 and onward LRI    UMR 8623 CNRS/Univ Paris Sud XI
///
///          Distributed under the Boost Software License, Version 1.0
///                 See accompanying file LICENSE.txt or copy at
///                     http://www.boost.org/LICENSE_1_0.txt
//////////////////////////////////////////////////////////////////////////////
/*!
 * \file
**/
#ifndef BOOST_SIMD_TOOLBOX_ARITHMETIC_FUNCTIONS_ISQRT_HPP_INCLUDED
#define BOOST_SIMD_TOOLBOX_ARITHMETIC_FUNCTIONS_ISQRT_HPP_INCLUDED
#include <boost/simd/include/simd.hpp>
#include <boost/dispatch/include/functor.hpp>

/*!
 * \ingroup boost_simd_arithmetic
 * \defgroup boost_simd_arithmetic_isqrt isqrt
 *
 * \par Description
 * return the integer truncation of the squre root of  the entry.
 * The value returned in the same type as the entry.
 * /par
 * Negative entry leads to zero for integer entry and nan for floating one.
 *
 * \par Header file
 * 
 * \code
 * #include <nt2/include/functions/isqrt.hpp>
 * \endcode
 * 
 * 
 * \synopsis
 *
 * \code
 * namespace boost::simd
 * {
 *   template <class A0>
 *     meta::call<tag::isqrt_(A0)>::type
 *     isqrt(const A0 & a0);
 * }
 * \endcode
 *
 * \param a0 the unique parameter of isqrt
 * 
 * \return an integer value
 *  
 * \par Notes
 * In SIMD mode, this function acts elementwise on the inputs vectors elements
 * \par
 *  
**/

namespace boost { namespace simd { namespace tag
  {         
    /*!
     * \brief Define the tag isqrt_ of functor isqrt 
     *        in namespace boost::simd::tag for toolbox boost.simd.arithmetic
    **/
    struct isqrt_ : ext::elementwise_<isqrt_> { typedef ext::elementwise_<isqrt_> parent; };
  }
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION(tag::isqrt_, isqrt, 1)
} }

#endif

// modified by jt the 25/12/2010
