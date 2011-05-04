//==============================================================================
//         Copyright 2003 & onward LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 & onward LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#ifndef NT2_SDK_ERROR_ERROR_HPP_INCLUDED
#define NT2_SDK_ERROR_ERROR_HPP_INCLUDED

/*!
 * \file
 * Implements NT2 exception handling system
 */

/*!
 * \defgroup error NT2 Error handling
 * \ingroup sdk
 * This module gathers macros, classes and functions to perform error handling
 * in NT2, beign at runtime or at compile-time.
 */

#include <nt2/sdk/error/warning.hpp>
#include <nt2/sdk/error/details/error.hpp>

////////////////////////////////////////////////////////////////////////////////
// No exception means no error unless they got requalified
////////////////////////////////////////////////////////////////////////////////
#include <boost/config.hpp>
#if defined(BOOST_NO_EXCEPTIONS)
NT2_WARNING(Exceptions globally disabled)
#define NT2_DISABLE_ERROR
#endif

////////////////////////////////////////////////////////////////////////////////
// Verbose report
////////////////////////////////////////////////////////////////////////////////
#if defined( NT2_VERBOSE )
  #if defined(NT2_CUSTOM_ERROR)
  NT2_WARNING(Using user-defined exception handler)
  #elif defined(NT2_DISABLE_ERROR)
  NT2_WARNING(Exceptions disabled)
  #endif
#endif

////////////////////////////////////////////////////////////////////////////////
// Enabled errors
////////////////////////////////////////////////////////////////////////////////
#if !defined(NT2_DISABLE_ERROR)
#include <nt2/sdk/error/details/exception.hpp>
#define NT2_THROW(EXP) BOOST_THROW_EXCEPTION( (EXP) ) \
/**/

////////////////////////////////////////////////////////////////////////////////
// Disabled errors
////////////////////////////////////////////////////////////////////////////////
#else
#define NT2_THROW(EXP)
#endif

#endif
