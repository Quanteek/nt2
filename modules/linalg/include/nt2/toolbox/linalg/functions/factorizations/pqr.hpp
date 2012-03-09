/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_FACTORIZATIONS_PQR_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_FACTORIZATIONS_PQR_HPP_INCLUDED

#include <nt2/toolbox/linalg/details/lapack/workspace.hpp>
#include <nt2/toolbox/linalg/details/utility/lapack_options.hpp>
#include <nt2/include/functions/pqr.hpp>
#include <nt2/include/functions/of_size.hpp>
#include <nt2/include/functions/min.hpp>
#include <nt2/include/functions/zeros.hpp>
#include <nt2/include/constants/one.hpp>
#include <nt2/include/functions/isempty.hpp>
#include <nt2/include/functions/leading_size.hpp>
#include <nt2/toolbox/linalg/details/lapack/pqr.hpp>
#include <nt2/toolbox/linalg/details/lapack/gqr.hpp>
#include <nt2/toolbox/linalg/details/lapack/mqr.hpp>
#include <nt2/toolbox/linalg/details/lapack/trs.hpp>

#include <nt2/table.hpp>
//#include <iostream>
//#include <nt2/include/functions/expand.hpp>
//#include <nt2/include/functions/diag.hpp>
//#include <nt2/include/functions/prod.hpp>
//#include <nt2/include/functions/triu.hpp>
//#include <nt2/include/functions/range.hpp>


//==============================================================================
// pqr call result forward declaration
//==============================================================================
namespace nt2
{
  template<class A> struct pqr_return;

} 

namespace nt2 { namespace ext
{
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::pqr_, tag::cpu_, 
                              (A)(S0),
                              ((expr_< table_<unspecified_<A>,S0>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                              )
  {
    typedef nt2::pqr_return<A> result_type; 
    BOOST_DISPATCH_FORCE_INLINE result_type operator()(const A& a) const
    {
      return nt2::pqr_return<A>(a);
    }
  };  
} }

namespace nt2 
{
  //============================================================================
  // pqr actual functor : precompute
  //============================================================================
  template<class A > struct pqr_return
  {
    typedef typename A::value_type                   type_t;
    typedef typename A::index_type                  index_t; 
    typedef typename meta::as_real<type_t>::type     base_t; 
    typedef nt2::table<type_t, nt2::matlab_index_>   ftab_t;
    typedef nt2::table<type_t, nt2::matlab_index_>  fbtab_t;
    typedef nt2::table<long int,nt2::matlab_index_>  iftab_t;
    typedef nt2::table<type_t, index_t>               tab_t;
    typedef nt2::table<type_t, index_t>              btab_t;
    typedef nt2::table<int32_t,index_t>              itab_t;
    typedef nt2::details::lapack_options          options_t; 
    typedef nt2::details::workspace<type_t>     workspace_t; 
    template < class XPR > pqr_return(const XPR& a_,
                                      workspace_t & w_):
      a(a_),
      m(size(a, 1)),
      n(size(a, 2)),
      k(std::min(m, n)), 
      jpvt(nt2::of_size(n, 1)), 
      lda(leading_size(a)), 
      tau(nt2::of_size(k, 1)), 
      mopts(), 
      w(w_)
    {
      init(); 
    }
    template < class XPR > pqr_return(const XPR & a_):
      a(a_),
      m(size(a, 1)),
      n(size(a, 2)),
      k(std::min(m, n)), 
      jpvt(nt2::of_size(n, 1)), 
      lda(leading_size(a)), 
      tau(nt2::of_size(k, 1)), 
      mopts(), 
      w(inw)
    {
      init(); 
    }
    ~pqr_return(){}
    // /////////////////////////////////////////////////////////////////////////////
    // accessors
    // /////////////////////////////////////////////////////////////////////////////
    tab_t       getq ()
    {
      if(is_empty(q))
        {
          q = a; //          q = expand(a, nt2::of_size(m, m));
          nt2::details::gqr(&m, &m, &k, q.raw(), &lda, tau.raw(), &info);
        }
      return q;
    }
    tab_t        getr()const
    {
      return triu(a);
    }
    tab_t       getp()
    {
      if(is_empty(p))
        {
          p = nt2::zeros(nt2::numel(jpvt), nt2::meta::as_<type_t>()); 
          for(unsigned int i=1; i <= nt2::size(p, 1) ; ++i){
            p(jpvt(i), i) = One<type_t>(); 
            }
          }
        return p; 
      }
      
    tab_t       getjp()         { return jpvt; }
    long int    getinfo() const { return info; }
    
    // /////////////////////////////////////////////////////////////////////////////
    size_t     rank(base_t epsi = nt2::Eps<base_t>())const
    {
      return 0; //globalSum(abs(diag(getr())) > nt2::max(n, m)*epsi*globalMax(abs(diag(getr()))) );
    }
    base_t absdet()const{
      //        mc_t::SquareTest(__FILE__, __LINE__, a);
      return 0; //prod(abs(diag(getr())));
    }
    
    // /////////////////////////////////////////////////////////////////////////////
    // resolutions
    // /////////////////////////////////////////////////////////////////////////////      
    template < class XPR > tab_t solve(const XPR & b, base_t epsi = nt2::Eps<base_t>(),
                                       const options_t & opts =  options_t()){
      ftab_t bb = b;
      
      const char side = 'L';
      const char tr = (opts.getTRANSA()) ? 'N' : !is_real(type_t(1))? 'C':'T';
      //        char tr = !isreal(type_t(1))? 'C':'T';
      const long int M = nt2::size(b, 1), N = nt2::size(b, 2); 
      const long int ldbb = nt2::leading_size(bb);
      nt2::details::mqr(&side, &tr, &M, &N, &k, a.raw(), &lda, tau.raw(), bb.raw(), &ldbb, &info, w);
      const long int nrhs = size(bb, 2); 
      const char uplo =  'U', d = 'N';
      const char tr1 = (opts.getTRANSA()) ?  !is_real(type_t(1))? 'C':'T' : 'N';
      //        tr =  'N'; 
      const long int rk = rank(epsi); 
      nt2::details::trtrs(&uplo, &tr1, &d, &rk, &nrhs, a.raw(), &lda, bb.raw(), &ldbb, &info, w);
      return permute(bb);
    }
      
    private :
      inline void init() 
      {
        nt2::details::geqp3(&m, &n, a.raw(), &lda, jpvt.raw(), tau.raw(), &info);
        //        mc_t::LapackTest(__FILE__, __LINE__, "geqp3", a, info); 
      }
      
    inline ftab_t permute(const ftab_t& bb) const {
      ftab_t res(nt2::of_size(nt2::numel(jpvt), nt2::size(bb, 2)));
      const size_t m =  nt2::min(size(bb, 1), numel(jpvt)); 
      //   res(jpvt(Range(1, m)), _) = bb; 
      return res; 
    }
    
    //      const char       jobu, jobvt;
    ftab_t                  a;
    ftab_t                  q;
    ftab_t                  p; 
    const long int    m, n, k;
    iftab_t              jpvt;
    const long int        lda; 
    ftab_t                tau; 
    long int             info;
    options_t           mopts; 
    workspace_t           inw; 
    workspace_t            &w;
  };
  
} 
#endif
