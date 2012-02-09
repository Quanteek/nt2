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
#ifndef NT2_TOOLBOX_ALGEBRA_FUNCTIONS_GEMM_HPP_INCLUDED
#define NT2_TOOLBOX_ALGEBRA_FUNCTIONS_GEMM_HPP_INCLUDED

#include <nt2/include/functor.hpp>

/*!
 * \ingroup algebra
 * \defgroup algebra_gemm gemm
 *
 * \par Description
 * Matrix multiplication
 *
 * \par Header file
 * 
 * \code
 * #include <nt2/include/functions/gemm.hpp>
 * \endcode
 * 
 * 
 * \synopsis
 *
 * \code
 * namespace nt2
 * {
 *   template <class A0,class A1>
 *     meta::call<tag::gemm_(A0,A1)>::type
 *     gemm(const A0 & a0, const A1 & a1);
 * }
 * \endcode
 *
 * \param a0 first parameter of the matrix product
 * 
 * \param a1 second parameter of the matrix product
 *
 * \return a value of the same type as the parameter
 *  
 * \par Notes
 * Call the dedicated blas routine available on the target.
 * \par
 *  
**/

namespace nt2 { namespace tag
  {         
    /*!
     * \brief Define the tag gemm_ of functor acos 
     *        in namespace nt2::tag for toolbox algebra
    **/
    struct gemm_ : ext::unspecified_<gemm_> { typedef ext::unspecified_<gemm_> parent; };
  }
  NT2_FUNCTION_IMPLEMENTATION_SELF(tag::gemm_, gemm, 3)
}

#endif
