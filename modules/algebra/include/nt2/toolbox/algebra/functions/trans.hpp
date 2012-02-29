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
#ifndef NT2_TOOLBOX_ALGEBRA_FUNCTIONS_TRANS_HPP_INCLUDED
#define NT2_TOOLBOX_ALGEBRA_FUNCTIONS_TRANS_HPP_INCLUDED

#include <nt2/include/functor.hpp>

/*!
 * \ingroup algebra
 * \defgroup algebra_trans trans
 *
 * \par Description
 * Matrix transposition
 *
 * \par Header file
 * 
 * \code
 * #include <nt2/include/functions/trans.hpp>
 * \endcode
 * 
 * 
 * \synopsis
 *     return transpose matrix with same settings but transposed dims and values
 *     if ta =  trans(a) ta(j, i) =  a(i, j) for every pair of correct indices
 *     for a.
 *     
 * \code
 * namespace nt2
 * {
 *   template <class A>
 *   auto trans(const &A a); 
 * }
 * \endcode
 *
 * \param a the matrix A
 *
 * \return  the transposed matrix
 *
 * \par
 *  
**/

namespace nt2 { namespace tag
  {         
    /*!
     * \brief Define the tag trans_ of functor trans
     *        in namespace nt2::tag for toolbox algebra
    **/
    struct trans_ : ext::unspecified_<trans_> { typedef ext::unspecified_<trans_> parent; };
  }

  template<class A>
  BOOST_FORCEINLINE void trans(A const& a)
  {
     typename nt2::make_functor<tag::trans_, A>::type()(a);
  }

}

#endif
