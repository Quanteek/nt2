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
 *     to the status chars transa and transb
 *     'N' meaning original
 *     'T' meaning transpose
 *     'H' meaning conjugate transpose
 * \code
 * namespace nt2
 * {
 *   template <class A0,class A1,class A2,class A3,class A4>
 *   void mm(const gem_status& gs,
 *             A0& c, A1 const& a, A2 & b,
 *             A3 const& alpha = 1, A4 const& beta = 0);
 *
 *   template <char transa,char transb,class A0,class A1,class A2,class A3,class A4>
 *   void mm(A0& c, A1 const& a, A2 & b,
 *             A3 const& alpha = 1, A4 const& beta = 0);
 * }
 * \endcode
 *
 * In template case:
 *
 * \param transa the transpose or no-transpose status of A: 'N' or 'T' or'H'
 *
 * \param transb the transpose or no-transpose status of B: 'N' or 'T' or'H'
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
     * \brief Define the tag mm_ of functor mm
     *        in namespace nt2::tag for toolbox algebra
    **/
    struct mm_ : ext::unspecified_<mm_> { typedef ext::unspecified_<mm_> parent; };
  }

  struct blas_types { //So this is not needed and is extracted from table layout ?
    static const std::size_t none                 = 0;
    static const std::size_t general              = 1;
    static const std::size_t band                 = 2;
    static const std::size_t packed               = 4;
    static const std::size_t symetric             = 8;
    static const std::size_t hermitian            = 16;
    static const std::size_t triangular           = 32; 
    static const std::size_t symetric_band        = symetric+band;
    static const std::size_t symetric_packed      = symetric+packed;
    static const std::size_t hermitian_band       = hermitian+band;
    static const std::size_t hermitian_packed     = hermitian+packed;
    static const std::size_t triangular_band      = triangular+band;
    static const std::size_t triangular_packed    = triangular+packed;
  };   

  template<std::size_t TYPE = blas_types::general,
           char UPLO   = 'L',
           char TRANSA = 'N',
           char TRANSB = 'N',
           char DIAG   = 'N',
           char SIDE   = 'L'>
  struct blas_status
  {
    static const int  type   = TYPE; 
    static const char uplo   = UPLO;
    static const char transa = TRANSA;
    static const char transb = TRANSB; 
    static const char diag   = DIAG;
    static const char side   = SIDE; 
  };


  template<class STATUS, class A0, class A1, class A2, class A3, class A4>
  //call with status
  BOOST_FORCEINLINE void
  mm(STATUS const& st, A0 const& a0, A1 const& a1, A2& a2, A3 const& a3, A4 const& a4)
  {
    return typename nt2::make_functor<tag::mm_, A0>::type()(st, a0, a1, a2, a3, a4);
  }
  template < class STATUS,  class A0,  class A1,  class A2,  class A3>
  BOOST_FORCEINLINE void mm(STATUS const& st, A0 const& a0, A1 const& a1, A2& a2, A3 const& a3)
  {
    typedef typename A0::value_type value_type; 
    mm(st, a0, a1, a2, a3, Zero<value_type>()); 
  }
  
  template < class STATUS,  class A0,  class A1,  class A2>
  BOOST_FORCEINLINE void mm(STATUS const& st, A0 const& a0, A1 const& a1, A2& a2)
  {
    typedef typename A0::value_type value_type; 
    mm(st, a0, a1, a2, One<value_type>(), Zero<value_type>()); 
  }

//   //templated on status
//   template < char transa,  char transb,  class A0,  class A1,  class A2,  class A3,  class A4>
//   BOOST_FORCEINLINE void mm(A0 const& a0, A1 const& a1, A2& a2, A3 const& a3, A4 const& a4)
//   {
//     typedef typename A0::value_type value_type; 
//     mm(mm_status<transa, transb>(), a0, a1, a2, a3, a4); 
//   }

//   template < int type, char transa,  char transb,  class A0,  class A1,  class A2,  class A3>
//   BOOST_FORCEINLINE void mm(A0 const& a0, A1 const& a1, A2& a2, A3 const& a3)
//   {
//     typedef typename A0::value_type value_type; 
//     mm(mm_status<transa, transb>(), a0, a1, a2, a3, Zero<value_type>()); 
//   }
  
//   template < int type,  char transa,  char transb, class A0,  class A1,  class A2>
//   BOOST_FORCEINLINE void mm(A0 const& a0, A1 const& a1, A2& a2)
//   {
//     typedef typename A0::value_type value_type; 
//     mm(mm_status<transa, transb>(), a0, a1, a2, One<value_type>(), Zero<value_type>()); 
//   }  
}

#endif
