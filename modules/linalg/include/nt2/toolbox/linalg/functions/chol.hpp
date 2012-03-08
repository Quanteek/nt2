/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2009-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_CHOL_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_CHOL_HPP_INCLUDED

#include <nt2/include/functor.hpp>

/*!
 * \ingroup algebra
 * \defgroup algebra_chol chol
 *
 * \par Description
 * Elementary Least square
 *
 * \par Header file
 * 
 * \code
 * #include <nt2/include/functions/chol.hpp>
 * \endcode
 * 
 * 
 * \synopsis
 * ggchol_itf.hh
    (excerpt adapted from xggchol.f file commentaries)
    
    DATA TYPE can mean float, double, std::complex<float>, std::complex<double>
    
    BASE TYPE can mean respectively float, double, float, double
    
    In some cases only two of these types types are available
    the two real or the two std::complex ones.
    CAPITALIZED PARAMETERS are FORTRAN parameters who are not used directly
    in the C++ calls, but through the workspace parameter,
    their use is transparent for the caller (see lapackworkspace.hh)

    *
    **  purpose
    **  =======
    **
    **  xggchol solves the linear equality-constrained least squares (chol)
    **  problem:
    **
    **          minimize || c - a*x ||_2   subject to   b*x = d
    **
    **  where a is an m-by-n matrix, b is a p-by-n matrix, c is a given
    **  m-vector, and d is a given p-vector. it is assumed that
    **  p <= n <= m+p, and
    **
    **           rank(b) = p and  rank( ( a ) ) = n.
    **                                ( ( b ) )
    **
    **  these conditions ensure that the chol problem has a unique solution,
    **  which is obtained using a grq factorization of the matrices b and a.
    **
    **  arguments (the C++ arguments are prefixed with #
    **  =========
    **
    **  m       (input) long int
    **          the number of rows of the matrix a.  m >= 0.
    **
    **  n       (input) long int
    **          the number of columns of the matrices a and b. n >= 0.
    **
    **  p       (input) long int
    **          the number of rows of the matrix b. 0 <= p <= n <= m+p.
    **
    **  #a       (input/output) DATA TYPE array, dimension (lda,n)
    **          on entry, the m-by-n matrix a.
    **          on exit, a is destroyed.
    **
    **  lda     (input) long int
    **          the leading dimension of the array a. lda >= max(1,m).
    **
    **  #b       (input/output) DATA TYPE array, dimension (ldb,n)
    **          on entry, the p-by-n matrix b.
    **          on exit, b is destroyed.
    **
    **  ldb     (input) long int
    **          the leading dimension of the array b. ldb >= max(1,p).
    **
    **  #c       (input/output) DATA TYPE array, dimension (m)
    **          on entry, c contains the right hand side vector for the
    **          least squares part of the chol problem.
    **          on exit, the residual sum of squares for the solution
    **          is given by the sum of squares of elements n-p+1 to m of
    **          vector c.
    **
    **  #d       (input/output) DATA TYPE array, dimension (p)
    **          on entry, d contains the right hand side vector for the
    **          constrained equation.
    **          on exit, d is destroyed.
    **
    **  #x       (output) DATA TYPE array, dimension (n)
    **          on exit, x is the solution of the chol problem.
    **
    **
    **
    **  info    (output) long int
    **          = 0:  successful exit.
    **          < 0:  if info = -i, the i-th argument had an illegal value.
    **

 *
 * \param a the matrix a on entry, destroyed on exit
 *
 * \param b the matrix b on entry, destroyed on exit
 * 
 * \param c the vector c on entry, residuals on exit
 *
 * \param d the vector d on entry, destroyed on exit
 *  
 * \param x  the solution of the chol problem on exit
 *
 * \par Notes
 *   Call the dedicated lapack routines available on the target.
 * \par
 *  
**/
//==============================================================================
// chol actual class forward declaration
//==============================================================================
namespace nt2 
{
  template<class A> struct chol_return;
} 

namespace nt2 { namespace tag
  {         
    /*!
     * \brief Define the tag chol_ of functor chol
     *        in namespace nt2::tag for toolbox algebra
    **/
    struct chol_ : ext::unspecified_<chol_> { typedef ext::unspecified_<chol_> parent; };
  }

  template<class A>
  BOOST_FORCEINLINE nt2::chol_return<A> chol(A &a)
  {
    return typename nt2::make_functor<tag::chol_, A>::type()(a);
  }

}

#endif

