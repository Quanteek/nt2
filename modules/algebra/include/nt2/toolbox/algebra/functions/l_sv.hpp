/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_ALGEBRA_FUNCTIONS_L_SV_HPP_INCLUDED
#define NT2_TOOLBOX_ALGEBRA_FUNCTIONS_L_SV_HPP_INCLUDED
#include <nt2/toolbox/algebra/lapack/lapackworkspace.hpp>
#include <nt2/toolbox/algebra/lapack/l_sv_Cpp.hpp>

////////////////////////////////////////////////////////////////////////////////
/*! \file l_sv.hpp
    **  sv computes the solution to a DATA TYPE system of linear equations
    **     a * x = b,
    **  where a is an n-by-n matrix and x and b are n-by-nrhs matrices.
    **
    **  the lu decomposition with partial pivoting and row interchanges is
    **  used to factor a as
    **     a = p * l * u,
    **  where p is a permutation matrix, l is unit lower triangular, and u is
    **  upper triangular.  the factored form of a is then used to solve the
    **  system of equations a * x = b.
**/

namespace nt2
{
  namespace ext
  {
    
    NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::l_sv_, tag::cpu_, 
                                (A)(SA)(B)(SB)(IPIV)(SIPIV)(STATUS), 
                                (unspecified_ < STATUS > )
                                ((table_< floating_<A>, SA > ))
                                ((table_< floating_<B>, SB > ))
                                ((table_< integer_<IPIV>, SIPIV > ))
                                )
    {
      typedef bool                                  result_type;
      typedef typename A::value_type                 value_type; 
      typedef lpp::Workspace<value_type>                    Wsp; 
      BOOST_FORCEINLINE result_type operator()( STATUS const&, A& a, B& b)
      {
        BOOST_ASSERT_MSG(is_square(a), "matrix a must be square in calling lapack sv routines"); 
        IPIV ipiv(of_size(size(a, 1))); 
        Wsp w; 
        return l_sv_call<A::shape>::l_sv_call( STATUS(), a, b, ipiv, w); 
      }
      BOOST_FORCEINLINE result_type operator()( STATUS const&, A& a, B& b, IPIV& ipiv)
      {
        BOOST_ASSERT_MSG(is_square(a), "matrix a must be square in calling lapack sv routines"); 
        lpp::Wsp w; 
        return l_sv_call<A::shape>::l_sv_call( STATUS(), a, b, ipiv, w); 
      }
      BOOST_FORCEINLINE result_type operator()( STATUS const&, A& a, B& b, IPIV&ipiv, Wsp& w)
      {
        BOOST_ASSERT_MSG(is_square(a), "matrix a must be square in calling lapack sv routines"); 
        return l_sv_call<A::form>::l_sv_call( STATUS(), a, b, ipiv, w); 
      }
      
    private:
      
      template < std::size_t FORM, class Dummy = void> struct blas_call
      {
        static inline result_type l_sv_call( STATUS const&, A&, B&, IPIV&, Wsp&)
        {
          BOOST_ASSERT_MSG(false, "The current matrix form is not supported by the lapack call");
          return false; 
        }
      };
      template < class Dummy = void> struct blas_call< l_type::general, Dummy>
      {
        static inline result_type l_sv_call( STATUS const&, A& a, B& b, IPIV& ipiv, Wsp&)
        {
          //           mcheck::SquareTest(__FILE__, __LINE__, A);
          //           mcheck::HeightsTest(__FILE__, __LINE__, A, X);
          long int info;
          const long int n     = nt2::extent(a)[0];
          const long int nrhs  = nt2::extent(x)[1];
          const long int lda   = nt2::details::padding(a);
          const long int ldb   = nt2::details::padding(b);
          gesv (&n, &nrhs, a.begin(), &lda, ipiv.begin(), x.begin(), &ldx, &info);
          //          mcheck::LapackTest(__FILE__, __LINE__, "gesv", A, info); 
          return (info == 0);
        }
      };

    }; 
}
#endif
