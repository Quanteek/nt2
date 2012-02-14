//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II       
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI       
//                                                                              
//          Distributed under the Boost Software License, Version 1.0.          
//                 See accompanying file LICENSE.txt or copy at                 
//                     http://www.boost.org/LICENSE_1_0.txt                     
//==============================================================================
/*!
 * \file
**/
#ifndef NT2_TOOLBOX_ALGEBRA_FUNCTIONS_GEMV_HPP_INCLUDED
#define NT2_TOOLBOX_ALGEBRA_FUNCTIONS_GEMV_HPP_INCLUDED

#include <nt2/include/functor.hpp>

/*!
 * \ingroup algebra
 * \defgroup algebra_gemv gemv
 *
 * \par Description
 * Matrix multiplication a0 <- alpha*a1*a2 + beta*a0
 * with 3 parameters alpha is 1 and beta 0
 * alpha and beta are scalars, a0, a1, a2 tables
 * of the same base type among float, double and
 * the associated complexes
 * If one parameter is complex all must be complex
 *
 * \par Header file
 * 
 * \code
 * #include <nt2/include/functions/gemv.hpp>
 * \endcode
 * 
 * 
 * \synopsis
 *
 * \code
 * namespace nt2
 * {
 *   template <class A0,class A1,class A2>
 *     meta::call<tag::gemv_(A0,A1,A2)>::type
 *     gemv(A0& a0, const A1 & a1, const A2 & a2);
 *   template <class A0,class A1,class A2,class A3>
 *     meta::call<tag::gemv_(A0,A1,A2,A3,A4)>::type
 *     gemv(A0& a0, const A1 & a1, const A2 & a2,
 *                  const A3 & alpha, const A4 & beta);
 * }
 * \endcode
 *
 * \param a0 first parameter of the matrix product in which the result is returned
 *
 * \param a1 second parameter of the matrix product
 * 
 * \param a2 third parameter of the matrix product  
 * 
 * \param alpha is optional and default to one
 *
 * \param beta is optional and default to zero
 *  
 * \par Notes
 * Call the dedicated blas routine available on the target.
 * \par
 *  
**/

namespace nt2 { namespace tag
  {         
    /*!
     * \brief Define the tag gemv_ of functor acos 
     *        in namespace nt2::tag for toolbox algebra
    **/
    struct general_gemv_ : ext::unspecified_<general_gemv_> { typedef ext::unspecified_<general_gemv_> parent; };
  }
  NT2_FUNCTION_IMPLEMENTATION_SELF(tag::general_gemv_, general_gemv, 6)

}

#endif
