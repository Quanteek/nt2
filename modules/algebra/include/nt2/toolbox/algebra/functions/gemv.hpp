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
 * Matrix vector multiplication
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
 *   template <class A0,class A1>
 *     meta::call<tag::gemv_(A0,A1,A2)>::type
 *     gemv(A0& a0, const A1 & a1, const A2 & a2);
 * }
 * \endcode
 *
 * \param a0 first parameter of the matrix vector product in which the result is returned
 * 
 * \param a1 second parameter of the matrix product, first matrix
 *
 * \param a2 second parameter of the matrix product, second matrix which is transposed
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
    struct gemv_ : ext::unspecified_<gemv_> { typedef ext::unspecified_<gemv_> parent; };
  }

  NT2_FUNCTION_IMPLEMENTATION_SELF(tag::gemv_, gemv, 3)
}

#endif
