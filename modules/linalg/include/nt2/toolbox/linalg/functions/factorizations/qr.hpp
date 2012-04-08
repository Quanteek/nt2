/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_FACTORIZATIONS_QR_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_FACTORIZATIONS_QR_HPP_INCLUDED

#include <nt2/toolbox/linalg/details/lapack/workspace.hpp>
#include <nt2/include/functions/qr.hpp>
#include <nt2/include/functions/of_size.hpp>
#include <nt2/include/functions/min.hpp>
#include <nt2/include/functions/isempty.hpp>
#include <nt2/include/functions/leading_size.hpp>
#include <nt2/toolbox/linalg/details/lapack/qr.hpp>
#include <nt2/toolbox/linalg/details/lapack/qrf.hpp>
#include <nt2/include/functions/triu.hpp>

#include <nt2/table.hpp>
//#include <iostream>
//#include <nt2/include/functions/expand.hpp>
//#include <nt2/include/functions/diag.hpp>
//#include <nt2/include/functions/prod.hpp>
//#include <nt2/include/functions/range.hpp>


//==============================================================================
// qr call result forward declaration
//==============================================================================
namespace nt2
{
  template<class A> struct qr_return;

} 

namespace nt2 { namespace ext
{
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::qr_, tag::cpu_, 
                              (A)(S0),
                              ((expr_< table_<unspecified_<A>,S0>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                              )
  {
    typedef nt2::qr_return<A> result_type; 
    BOOST_DISPATCH_FORCE_INLINE result_type operator()(const A& a) const
    {
      return nt2::qr_return<A>(a);
    }
  };  
} }

namespace nt2 
{
  //============================================================================
  // qr actual functor : precompute
  //============================================================================
  template<class A > struct qr_return
  {
    typedef typename A::value_type                   type_t;
    typedef typename A::index_type                  index_t; 
    typedef typename meta::as_real<type_t>::type     base_t; 
    typedef nt2::table<type_t, nt2::matlab_index_>   ftab_t;
    typedef nt2::table<type_t, nt2::matlab_index_>  fbtab_t;
    typedef nt2::table<type_t, index_t>               tab_t;
    typedef nt2::table<type_t, index_t>              btab_t;
    
    template < class XPR > qr_return(const XPR & a_):
      a(a_),
      m(size(a, 1)),
      n(size(a, 2)),
      k(std::min(m, n)), 
      lda(leading_size(a)), 
      tau(nt2::of_size(k, 1)), 
      w(inw)
    {
      init(); 
    }
    template < class XPR > qr_return(const XPR & a_,
                                     nt2::details::workspace < type_t > & w_):
      a(a_),
      m(size(a, 1)),
      n(size(a, 2)),
      k(std::min(m, n)), 
      lda(leading_size(a)), 
      tau(nt2::of_size(k, 1)), 
      w(w_)
    {
      init(); 
    }
    ~qr_return(){}
    // /////////////////////////////////////////////////////////////////////////////
    // accessors
    // /////////////////////////////////////////////////////////////////////////////
    tab_t      getq()
    {
      if (isempty(q)){
        long int nn = nt2::min(n, m); //length(a); 
        q = a; //         q = expand(a, nt2::of_size(m, nn));
        nt2::details::gqr(&m, &nn, &k, q.raw(), &m, tau.raw(), &info, w);
      }
      return q;
    }
    tab_t       getr()    const { return triu(a/*(Range(0, k-1), _)*/); }
    long int    getinfo() const { return info; }
    
    // /////////////////////////////////////////////////////////////////////////////
    //       size_t     rank(base_t epsi = eps<base_t>())
    //       {
    //         return globalSum( abs(diag(getr())) > std::max(n, m)*epsi*globalMax(abs(diag(getr()))) );
    //       }
    //       base_t absdet(){
    //         //        mc_t::SquareTest(__FILE__, __LINE__, a);
    //         return nt2::prod(nt2::abs(nt2::diag(nt2::getr())));
    //       }
    
    
    // /////////////////////////////////////////////////////////////////////////////
    // resolutions
    // /////////////////////////////////////////////////////////////////////////////
    
  private :
    inline void init()
    {
      nt2::details::geqrf(&m, &n, a.raw(), &lda, tau.raw(), &info, w);
      //        mc_t::LapackTest(__FILE__, __LINE__, "geqrf", a, info); 
    }
    //      const char     jobu, jobvt;
    tab_t                  a;
    tab_t                  q;
    const long int   m, n, k;
    const long int       lda; 
    tab_t                tau; 
    long int            info; 
    nt2::details::workspace < type_t > inw; 
    nt2::details::workspace < type_t >  &w;
  };
} 
#endif
