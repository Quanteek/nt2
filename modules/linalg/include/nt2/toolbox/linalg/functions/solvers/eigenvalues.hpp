/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_SOLVERS_EIGENVALUES_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_SOLVERS_EIGENVALUES_HPP_INCLUDED

#include <nt2/include/functions/eigenvalues.hpp>
#include <nt2/include/functions/numel.hpp>
#include <nt2/include/functions/leading_size.hpp>
#include <nt2/toolbox/linalg/details/lapack/ev.hpp>
#include <nt2/table.hpp>
#include <nt2/include/functions/max.hpp>
#include <nt2/include/functions/issquare.hpp>
#include <nt2/include/functions/ishermitian.hpp>
#include <nt2/include/functions/height.hpp>
#include <nt2/include/functions/width.hpp>


//==============================================================================
// eigenvalues actual functor forward declaration
//==============================================================================
namespace nt2
{
  template<class A, class V, class B, class C> struct eigenvalues_return;
} 

namespace nt2 { namespace ext
{
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::eigenvalues_, tag::cpu_, 
                              (A)(SA)(V)(SV)(B)(C),
                              ((expr_< table_<unspecified_<A>,SA>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                              ((expr_< table_<unspecified_<V>,SV>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                              (scalar_<integer_<B> >
                              (scalar_<integer_<C> >
                              )
  {
    typedef nt2::eigenvalues_return<A, B, C> result_type; 
    BOOST_DISPATCH_FORCE_INLINE result_type operator()(A& a, V& v, const B& b, const C& c) const
    {
      return nt2::eigenvalues_return<A>(a, v, b, c);
    }
  };  
} }

namespace nt2 
{
  //============================================================================
  // svd actual functor : precompute
  //============================================================================
  template<class A, class V> struct eigenvalues_return
  {
    typedef typename A::value_type                   type_t;
    typedef typename A::index_type                  index_t; 
    typedef typename meta::as_real<type_t>::type    btype_t; 
    typedef nt2::table<type_t,nt2::matlab_index_>    ftab_t;
    typedef nt2::table<btype_t,nt2::matlab_index_>  fbtab_t;
    typedef nt2::table<la_int,nt2::matlab_index_>   fitab_t;
    typedef nt2::table<type_t,index_t>                tab_t;
    typedef nt2::table<btype_t,index_t>              btab_t;
    typedef nt2::table<la_int,index_t>               itab_t;

    eigenvalues_return(A& a,  V& v, const char & jobz = 'n', const char& uplo = 'l')
    {
      if (lower(jobz) !=  'n' && lower(jobz) != 'l')
        {
          BOOST_ASSERT_MSG(nt2::ishermitian(a), "matrix a is not hermitian");
          uplo =  'l';
        }
      BOOST_ASSERT_MSG(nt2::issquare(a), "matrix a is not square");
      const la_int n = height(a);
      const la_int lda = nt2::leading_size(a);
      v.resize(nt2::of_size(height(a), 1));
      ev(&jobz, &uplo, &n, a.raw(), &lda, v.raw(), &info);
    }
    la_int get_info()    const { return info; }
  private:
    la_int  info;
  };
} 

#endif
