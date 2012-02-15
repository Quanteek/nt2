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
 *   template <class A0,class A1,class A1, class A2, class A3, class A4, class A5>
 *     meta::call<tag::gemv_(A5,A0,A1,A2,A3,A4)>::type
 *     gemv(A5 const& a5, A0 const& a0, A1 const& a1, A2 & a2, A3 const& a3, A4 const& a4);
 * }
 * \endcode
 *
 * \param a0 first parameter of the matrix-vector product, the input matrix
 *
 * \param a1 second parameter of the matrix-vector product, the input vector
 * 
 * \param a2 third parameter of the matrix-vector product in which the result is returned
 *  
 * \param a3 fourth parameter of the matrix-vector product, alpha parameter
 *
 * \param a4 fifth parameter of the matrix-vector product, beta parameter
 *
 * \param a5 specifies the transpose or none-transpose parameter with the class gemm_status
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

  template<class A0, class A1, class A2, class A3, class A4, class A5>
  BOOST_FORCEINLINE 
  typename nt2::meta::
  call<nt2::tag::gemv_( A5 const&, A0 const&, A1 const&
                      , A2&      , A3 const&
                      , A4 const&)
      >::type
  gemv(A5 const& a5, A0 const& a0, A1 const& a1, A2& a2, A3 const& a3, A4 const& a4)
  {
    return typename nt2::make_functor<tag::gemv_, A0>::type()(a5, a0, a1, a2, a3, a4);
  }

}

#endif
