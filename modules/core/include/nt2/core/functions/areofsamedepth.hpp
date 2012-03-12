//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#ifndef NT2_CORE_FUNCTIONS_AREOFSAMEDEPTH_HPP_INCLUDED
#define NT2_CORE_FUNCTIONS_AREOFSAMEDEPTH_HPP_INCLUDED

#include <nt2/include/functor.hpp>

/*!
 * \ingroup core
 * \defgroup core_areofsamedepth areofsamedepth
 *
 * \par Description
 * Returns true or false according a0 and a1 are of same depth
 *
 * \par Header file
 * 
 * \code
 * #include <nt2/include/functions/areofsamedepth.hpp>
 * \endcode
 * 
 * \par Alias 
 * \arg are_of_same_depth
 * 
 * \synopsis
 *
 * \code
 * namespace boost::simd
 * {
 *   template <class A0, class A1>
 *     bool areofsamedepth(const A0 & a0, const A1 & a1);
 * }
 * \endcode
 *
 * \param a0 the first parameter of areofsamedepth
 * \param a1 the second parameter of areofsamedepth
 * 
 * \return a bool value
 *  
**/

namespace nt2
{
  namespace tag
  {
    struct areofsamedepth_ : ext::unspecified_<areofsamedepth_>
    {
      typedef ext::unspecified_<areofsamedepth_> parent;
    };
  }

  NT2_FUNCTION_IMPLEMENTATION(nt2::tag::areofsamedepth_, areofsamedepth, 2)
  NT2_FUNCTION_IMPLEMENTATION(nt2::tag::areofsamedepth_, are_of_same_depth, 2)
}

#endif

