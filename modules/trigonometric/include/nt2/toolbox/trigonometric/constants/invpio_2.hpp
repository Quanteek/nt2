//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II       
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI       
//                                                                              
//          Distributed under the Boost Software License, Version 1.0.          
//                 See accompanying file LICENSE.txt or copy at                 
//                     http://www.boost.org/LICENSE_1_0.txt                     
//==============================================================================
#ifndef NT2_TOOLBOX_TRIGONOMETRIC_CONSTANTS_INVPIO_2_HPP_INCLUDED
#define NT2_TOOLBOX_TRIGONOMETRIC_CONSTANTS_INVPIO_2_HPP_INCLUDED
/*!
 * \file
**/
#include <nt2/sdk/constant/constant.hpp>
/*!
 * \ingroup boost_simd_constant
 * \defgroup boost_simd_constant_invpio_2 invpio_2 constant
 *
 * \par Description
 * Constant Invpio_2 : \f$\frac1\pi\f$.
 * 
 * \par Header file
 * 
 * \code
 * #include <nt2/include/functions/invpio_2.hpp>
 * \endcode
 * 
 * 
 * \synopsis
 *
 * \code
 * namespace boost::simd
 * {
 *   template <class T,class A0>
 *     meta::call<tag::_invpio_2_(A0)>::type
 *     Invpio_2();
 * }
 * \endcode
 *
 * 
 * \param T template parameter of Invpio_2
 * 
 * \return type T value
 *  
 *  
**/
namespace nt2
{
  NT2_MAKE_STD_CONSTANT(Invpio_2   , 0x3FE45F306DC9C883ll, 0x3f22f984, 0)//6.36619772367581382433e-01
}

#endif
