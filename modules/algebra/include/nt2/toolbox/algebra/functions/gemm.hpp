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
 *     matrix product encapsulasion of blas (dszc)gemm
 *     C <- alpha*A°*B°+beta*C
 *     A and B are matrices and A° and B° are (conjugate) transpose, or original according
 *     to the status letters transa and transb
 *     'N' meaning original
 *     'T' meaning transpose
 *     'H' meaning conjugate transpose
 * \code
 * namespace nt2
 * {
 *   template <class A0,class A1,class A2,class A3,class A4>
 *   void gemm(const gem_status& gs,
 *             A0& c, A1 const& a, A2 & b,
 *             A3 const& alpha = 1, A4 const& beta = 0);
 *
 *   template <char transa,char transb,class A0,class A1,class A2,class A3,class A4>
 *   void gemm(A0& c, A1 const& a, A2 & b,
 *             A3 const& alpha = 1, A4 const& beta = 0);
 * }
 * \endcode
 *
 * In template case:
 *
 * \param transa the transpose or no-transpose status of A
 *
 * \param transb the transpose or no-transpose status of B
 *
 * In call with status case: the first parameter gs is simply nt2::gem_status<transa,transb>()
 *
 * Other parameters are the same for both type of calls:
 *
 * \param c the matrix product result and C input
 *
 * \param a the matrix A
 * 
 * \param b the matrix B
 *  
 * \param alpha the alpha scalar parameter
 *
 * \param beta  the beta scalar parameter
 *
 * \par Notes
 *   Call the dedicated blas routines available on the target.
 * \par
 *  
**/

namespace nt2 { namespace tag
  {         
    /*!
     * \brief Define the tag gemm_ of functor gemm
     *        in namespace nt2::tag for toolbox algebra
    **/
    struct gemm_ : ext::unspecified_<gemm_> { typedef ext::unspecified_<gemm_> parent; };
  }

  template<class A0, class A1, class A2, class A3, class A4, class A5>
  BOOST_FORCEINLINE 
  typename nt2::meta::
  call<nt2::tag::gemm_( A5 const&, A0 const&, A1 const&
                      , A2&      , A3 const&
                      , A4 const&)
      >::type
  gemm(A5 const& a5, A0 const& a0, A1 const& a1, A2& a2, A3 const& a3, A4 const& a4)
  {
    return typename nt2::make_functor<tag::gemm_, A0>::type()(a5, a0, a1, a2, a3, a4);
  }
}

#endif
