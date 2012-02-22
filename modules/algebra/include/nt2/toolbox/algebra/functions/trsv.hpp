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
#ifndef NT2_TOOLBOX_ALGEBRA_FUNCTIONS_TRSV_HPP_INCLUDED
#define NT2_TOOLBOX_ALGEBRA_FUNCTIONS_TRSV_HPP_INCLUDED

#include <nt2/include/functor.hpp>

/*!
 * \ingroup algebra
 * \defgroup algebra_trsv trsv
 *
 * \par Description
 * Triangular Matrix vector multiplication
 *
 * \par Header file
 * 
 * \code
 * #include <nt2/include/functions/trsv.hpp>
 * \endcode
 * 
 * 
 * \synopsis
 *     matrix product encapsulasion of blas (dszc)trsv
 *    y <- inv(A°)x°
 *    A° is original, transposed or transcongugated A, according transa 'N', 'T' or 'H'
 *    x° is original or transposed x, according transx 'N', 'T' or 'H'
 *    if upload is 'U', A is supposed to be upper-triangular else lower-triangular
 *    if diag   is 'U', A is supposed to have ones on the diagonal
 *    The two last information are taken as ground truth an the algorithm acts as they
 *    were truly verified. However the non relevant values in A are nor accessed, nor modified.
 * 
 * \code
 * namespace nt2
 * {
 *   template <class A0,class A1,class A2>
 *   void gemm(const trsv_status& gs,
 *             A0& a, A1 const& x, A2 & y);
 *
 *   template <char transa,char transx, char upload, char diag,
 *             class A0,class A1,class A2>
 *   void gemm(A0& a, A1 const& x, A2 & y,
 *             A3 const& alpha = 1, A4 const& beta = 0);
 * }
 * \endcode
 *
 * In template case:
 *
 * \param transa the transpose or no-transpose status of A: 'N' or 'T' or 'C'
 *
 * \param transx the transpose or no-transpose status of x: 'N', never conjugated
 *
 * In call with status case: the first parameter gs is simply nt2::trsv_status<transa,transx>()
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
     * \brief Define the tag trsv_ of functor trsv 
     *        in namespace nt2::tag for toolbox algebra
    **/
    struct trsv_ : ext::unspecified_<trsv_> { typedef ext::unspecified_<trsv_> parent; };
  }

  template<char T0, char T1, char T2, char T3>
  struct trsv_status
  {
    static const char transa = T0;
    static const char transx = T1;
    static const char upload = T2;
    static const char diag   = T3; 
  };

  template<class A0, class A1, class A5>
  BOOST_FORCEINLINE void
  trsv(A5 const& a5, A0 const& a0, A1 & a1)
  {
    typename nt2::make_functor<tag::trsv_, A0>::type()(a5, a0, a1);
  }

  template < char transa, char transx, char upload, char diag,
             class A0, class A1, class A2>
  BOOST_FORCEINLINE void trsv(A0 const& a0, A1& a1)
  {
    trsv(trsv_status<transa, transx, upload, diag>(), a0, a1);
  }
}

#endif
