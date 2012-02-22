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
 *   template <class A0,class A1,class A2,class A3,class A4>
 *   void gemm(const gemv_status& gs,
 *             A0& a, A1 const& x, A2 & y,
 *             A3 const& alpha = 1, A4 const& beta = 0);
 *
 *   template <char transa,char transx,class A0,class A1,class A2,class A3,class A4>
 *   void gemm(A0& a, A1 const& x, A2 & y,
 *             A3 const& alpha = 1, A4 const& beta = 0);
 * }
 * \endcode
 *
 * In template case:
 *
 * \param transa the transpose or no-transpose status of A: 'N' or 'T' or 'H'
 *
 * \param transx the transpose or no-transpose status of x: 'N', never conjugated
 *
 * In call with status case: the first parameter gs is simply nt2::gemv_status<transa,transx>()
 *
 * Other parameters are the same for both type of calls:
 *
 * \param a the matrix A
 * 
 * \param b the vector x
 *  
 * \param y the column vector product result and input
 *
 * \param alpha the alpha scalar parameter
 *
 * \param beta  the beta scalar parameter
 *
 * \par Notes
 * Call the dedicated blas routine available on the target.
 * \par
 *  
**/

namespace nt2 { namespace tag
  {         
    /*!
     * \brief Define the tag gemv_ of functor gemv 
     *        in namespace nt2::tag for toolbox algebra
    **/
    struct gemv_ : ext::unspecified_<gemv_> { typedef ext::unspecified_<gemv_> parent; };
  }

  template<char T0, char T1>
  struct gemv_status
  {
    static const char tA = T0;
    static const char tx = T1;
  };

  template<class A0, class A1, class A2, class A3, class A4, class A5>
  BOOST_FORCEINLINE void
  gemv(A5 const& a5, A0 const& a0, A1 const& a1, A2& a2, A3 const& a3, A4 const& a4)
  {
    return typename nt2::make_functor<tag::gemv_, A0>::type()(a5, a0, a1, a2, a3, a4);
  }
 
  template < class A5, class A0, class A1, class A2, class A3>
  BOOST_FORCEINLINE void gemv(A5 const& a5, A0 const& a0, A1 const& a1, A2& a2, A3 const& a3)
  {
    typedef typename A0::value_type value_type; 
    gemv(a5, a0, a1, a2, a3, Zero<value_type>());
  }
  
  template < class A5, class A0, class A1, class A2>
  BOOST_FORCEINLINE void gemv(A5 const& a5, A0 const& a0, A1 const& a1, A2& a2)
  {
    typedef typename A0::value_type value_type; 
    gemv(a5, a0, a1, a2, One<value_type>(), Zero<value_type>());
  }
  template < char transa, char transx, class A0, class A1, class A2, class A3, class A4>
  BOOST_FORCEINLINE void gemv(A0 const& a0, A1 const& a1, A2& a2, A3 const& a3, A4 const& a4)
  {
    gemv(gemv_status<transa, transx>(), a0, a1, a2, a3, a4);
  }
  
  template < char transa, char transx, class A0, class A1, class A2, class A3>
  BOOST_FORCEINLINE void gemv(A0 const& a0, A1 const& a1, A2& a2, A3 const& a3)
  {
    typedef typename A0::value_type value_type; 
    gemv(gemv_status<transa, transx>(), a0, a1, a2, a3, Zero<value_type>());
  }
  
  template <char transa, char transx, class A0, class A1, class A2>
  BOOST_FORCEINLINE void gemv(A0 const& a0, A1 const& a1, A2& a2)
  {
    typedef typename A0::value_type value_type; 
    gemv(gemv_status<transa, transx>(), a0, a1, a2, One<value_type>(), Zero<value_type>());
  }
  
}

#endif
