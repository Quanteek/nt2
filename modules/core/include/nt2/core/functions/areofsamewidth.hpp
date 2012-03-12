//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#ifndef NT2_CORE_FUNCTIONS_AREOFSAMEWIDTH_HPP_INCLUDED
#define NT2_CORE_FUNCTIONS_AREOFSAMEWIDTH_HPP_INCLUDED

#include <nt2/include/functor.hpp>

/*!
 * \ingroup core
 * \defgroup core_areofsamewidth areofsamewidth
 *
 * \par Description
 * Returns true or false according a0 and a1 are of same width
 *
 * \par Header file
 * 
 * \code
 * #include <nt2/include/functions/areofsamewidth.hpp>
 * \endcode
 * 
 * \par Alias 
 * \arg are_of_same_width
 * 
 * \synopsis
 *
 * \code
 * namespace boost::simd
 * {
 *   template <class A0, class A1>
 *     bool areofsamewidth(const A0 & a0, const A1 & a1);
 * }
 * \endcode
 *
 * \param a0 the first parameter of areofsamewidth
 * \param a1 the second parameter of areofsamewidth
 * 
 * \return a bool value
 *  
**/

namespace nt2
{
  namespace tag
  {
    struct areofsamewidth_ : ext::unspecified_<areofsamewidth_>
    {
      typedef ext::unspecified_<areofsamewidth_> parent;
    };
  }

  NT2_FUNCTION_IMPLEMENTATION(nt2::tag::areofsamewidth_, areofsamewidth, 2)
  NT2_FUNCTION_IMPLEMENTATION(nt2::tag::areofsamewidth_, are_of_same_width, 2)
}

#endif

