/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2009-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_SV_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_SV_HPP_INCLUDED

#include <nt2/include/functor.hpp>

/*!
 * \ingroup algebra
 * \defgroup algebra_sv sv
 *
 * \par Description
 * Elementary Least square
 *
 * \par Header file
 * 
 * \code
 * #include <nt2/include/functions/sv.hpp>
 * \endcode
 * 
 * 
 * \synopsis
 *  sv computes the solution to a DATA TYPE system of linear equations
 *     a * x = b,
 *  where a is an n-by-n matrix and x and b are n-by-nrhs matrices.
 *
 *  the lu decomposition with partial pivoting and row interchanges is
 *  used to factor a as
 *     a = p * l * u,
 *  where p is a permutation matrix, l is unit lower triangular, and u is
 *  upper triangular.  the factored form of a is then used to solve the
 *  system of equations a * x = b.
 *
 * \param a the matrix a on entry, destroyed on exit
 *
 * \param b the matrix b on entry, result x on exit
 * 
 * \par Notes
 *   Call the dedicated lapack routines available on the target.
 * \par
 *  
**/

namespace nt2 { namespace tag
  {         
    /*!
     * \brief Define the tag sv_ of functor sv
     *        in namespace nt2::tag for toolbox algebra
    **/
    struct sv_ : ext::unspecified_<sv_> { typedef ext::unspecified_<sv_> parent; };
  }

  template<class A, class B>
  BOOST_FORCEINLINE long int sv(A &a, B& b)
  {
    return typename nt2::make_functor<tag::sv_, A>::type()(a, b);
  }
  template<class A, class B, class IPIV>
  BOOST_FORCEINLINE long int sv(A &a, B& b, IPIV& ipiv)
  {
     return typename nt2::make_functor<tag::sv_, A>::type()(a, b, ipiv);
  }

}

#endif

