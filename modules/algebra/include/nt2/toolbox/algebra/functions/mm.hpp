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
#ifndef NT2_TOOLBOX_ALGEBRA_FUNCTIONS_MM_HPP_INCLUDED
#define NT2_TOOLBOX_ALGEBRA_FUNCTIONS_MM_HPP_INCLUDED

#include <nt2/include/functor.hpp>

/*!
 * \ingroup algebra
 * \defgroup algebra_mm mm
 *
 * \par Description
 * Matrix multiplication
 *
 * \par Header file
 * 
 * \code
 * #include <nt2/include/functions/mm.hpp>
 * \endcode
 * 
 * 
 * \synopsis
 *     matrix product encapsulasion of blas (d, s, z, c)(ge, gb, sy, he)mm
 *     C <- alpha*A°*B°+beta*C
 *     A and B are matrices and A° and B° are (conjugate) transpose, or original according
 *     to the status chars TRANSA and TRANSB, SIDE, DIAG and UPLO
 *     For TRANSA and TRANSB (no effects for the 
 *       'N' meaning no transpose
 *       'T' meaning transpose
 *       'C' meaning conjugate transpose
 *     
 *     Note that 
 * \code
 * namespace nt2
 * {
 *   template <class A0,class A1,class A2,class A3,class A4>
 *   void mm(const blas_status& gs,
 *             A0& c, A1 const& a, const A2& b,
 *             A3 const& alpha = 1, A4 const& beta = 0);
 * }
 * \endcode
 *
 * the first parameter bs 
 *    nt2::blas_status<UPLO, TRANSA, TRANSB, DIAG, SIDE>()
 * where the five template parameters are chars defining the blas status.
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
     * \brief Define the tag mm_ of functor mm
     *        in namespace nt2::tag for toolbox algebra
    **/
    struct mm_ : ext::unspecified_<mm_> { typedef ext::unspecified_<mm_> parent; };
  }

  template<char UPLO   = 'L',
           char TRANSA = 'N',
           char TRANSB = 'N',
           char DIAG   = 'N',
           char SIDE   = 'L'>
  struct b_status
  {
    static const char uplo   = UPLO;
    static const char transa = TRANSA;
    static const char transb = TRANSB; 
    static const char diag   = DIAG;
    static const char side   = SIDE; 
  };


  template<class STATUS, class A, class B, class C, class Alpha, class Beta>
  BOOST_FORCEINLINE void mm(STATUS const& st, A const& a, B const& b,
                            C& c,
                            Alpha const& alpha, Beta const& beta)
  {
     typename nt2::make_functor<tag::mm_, A>::type()(st, a, b, c, alpha, beta);
  }

  template < class STATUS,  class A, class B, class C, class Alpha>
  BOOST_FORCEINLINE void mm(STATUS const& st, A const& a, B const& b,
                            C& c,
                            Alpha const& alpha)
  {
    typedef typename C::value_type value_type; 
    mm(st, a, b, c, alpha, Zero<value_type>()); 
  }
  
  template < class STATUS,  class A, class B, class C>
  BOOST_FORCEINLINE void mm(STATUS const& st, A const& a, B const& b,
                            C& c)
  {
    typedef typename C::value_type value_type; 
    mm(st, a, b, c, One<value_type>(), Zero<value_type>()); 
  }
  
}

#endif
