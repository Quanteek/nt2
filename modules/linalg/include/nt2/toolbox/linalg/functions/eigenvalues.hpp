/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2009-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_EIGENVALUES_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_EIGENVALUES_HPP_INCLUDED
#include <nt2/include/functor.hpp>

/*!
 * \ingroup algebra
 * \defgroup algebra_eigenvalues eigenvalues
 *
 * \par Description
 * Elementary Least square
 *
 * \par Header file
 * 
 * \code
 * #include <nt2/include/functions/eigenvalues>
 * \endcode
 * 
 * 
 * \synopsis
 * eigen values of a symetric matrix
 * if uplo is 'L' or 'U' not check is made about the fact a is really hermitian
 * else the symetry is checked.
 *
 * \param a the matrix a on entry, replaced on exit by eigen vectors if jobz is 'V'
 *                                 untouched if jobz is 'N'
 *
 * \param values undefined on entry, eigen values on exit
 * 
 * \param jobz  'V' or 'N',  if 'N' eigenvectors are not computed
 
 * 
 * \par Notes
 *   Call the dedicated lapack routines available on the target.
 * \par
 *  
**/
//==============================================================================
// eigenvalues actual class forward declaration
//==============================================================================
namespace nt2 
{
  template<class A, class V, class B, class C> struct eigenvalues_return;
} 

namespace nt2 { namespace tag
  {         
    /*!
     * \brief Define the tag eigenvalues_ of functor eigenvalues
     *        in namespace nt2::tag for toolbox algebra
    **/
    struct eigenvalues_ : ext::unspecified_<eigenvalues_> { typedef ext::unspecified_<eigenvalues_> parent; };
  }

  template<class A, class V>
  BOOST_FORCEINLINE nt2::eigenvalues_return<A, V> eigenvalues(A &a, V& values, const char & jobz = 'v', const char uplo = 'v')
  {
    return typename nt2::make_functor<tag::eigenvalues_, A>::type()(a, jobz, uplo);
  }

}

#endif

