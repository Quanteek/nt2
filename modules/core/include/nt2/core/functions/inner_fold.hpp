//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#ifndef NT2_CORE_FUNCTIONS_INNER_FOLD_HPP_INCLUDED
#define NT2_CORE_FUNCTIONS_INNER_FOLD_HPP_INCLUDED

/*!
 * \file
 * \brief Defines and implements the nt2::inner_fold function
 */

#include <nt2/include/functor.hpp>

namespace nt2
{
  namespace tag { struct inner_fold_ : boost::dispatch::tag::formal_ { typedef boost::dispatch::tag::formal_ parent; }; }

  //============================================================================
  /*!
   * Folds elements of \c a1 along inner dimension, possibly in parallel, and store
   * the result in \c a0.
   *
   * \param a0 Expression to store result in
   * \param a1 Expression to reduce
   * \param a2 Functor to initialize with
   * \param a3 Function to apply for reduction (binary + unary if SIMD is used)
   * \return nothing
   */
  //============================================================================
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION_TPL(tag::inner_fold_, inner_fold, (A0 const&)(A1 const&)(A2 const&)(A3 const&)(A4 const&), 5)
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION_TPL(tag::inner_fold_, inner_fold, (A0&)(A1 const&)(A2 const&)(A3 const&)(A4 const&), 5)
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION_TPL(tag::inner_fold_, inner_fold, (A0 const&)(A1&)(A2 const&)(A3 const&)(A4 const&), 5)
  BOOST_DISPATCH_FUNCTION_IMPLEMENTATION_TPL(tag::inner_fold_, inner_fold, (A0&)(A1&)(A2 const&)(A3 const&)(A4 const&), 5)

}



#endif
