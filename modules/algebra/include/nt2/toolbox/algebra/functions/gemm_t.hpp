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
#ifndef NT2_TOOLBOX_ALGEBRA_FUNCTIONS_GEMM_T_HPP_INCLUDED
#define NT2_TOOLBOX_ALGEBRA_FUNCTIONS_GEMM_T_HPP_INCLUDED

#include <nt2/include/functor.hpp>

/*!
 * \ingroup algebra
 * \defgroup algebra_gemm_t gemm_t
 *
 * \par Description
 * Matrix multiplication
 *
 * \par Header file
 * 
 * \code
 * #include <nt2/include/functions/gemm_t.hpp>
 * \endcode
 * 
 * 
 * \synopsis
 *
 * \code
 * namespace nt2
 * {
 *   template <class A0,class A1>
 *     meta::call<tag::gemm_t_(A0,A1,A2)>::type
 *     gemm_t(A0& a0, const A1 & a1, const A2 & a2);
 * }
 * \endcode
 *
 * \param a0 first parameter of the matrix product in which the result is returned
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
     * \brief Define the tag gemm_t_ of functor acos 
     *        in namespace nt2::tag for toolbox algebra
    **/
    struct gemm_t_ : ext::unspecified_<gemm_t_> { typedef ext::unspecified_<gemm_t_> parent; };
  }
  NT2_FUNCTION_IMPLEMENTATION_SELF(tag::gemm_t_, gemm_t, 3)
}

#endif
