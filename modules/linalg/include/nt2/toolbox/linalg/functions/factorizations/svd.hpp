/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_FACTORIZATIONS_SVD_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_FACTORIZATIONS_SVD_HPP_INCLUDED

#include <nt2/include/functions/svd.hpp>
#include <nt2/include/functions/numel.hpp>
#include <nt2/include/functions/leading_size.hpp>
#include <nt2/toolbox/linalg/details/utility/tags.hpp>
#include <nt2/toolbox/linalg/details/lapack/workspace.hpp>
#include <nt2/toolbox/linalg/details/utility/lower.hpp>
#include <nt2/toolbox/linalg/details/lapack/svd.hpp>
#include <nt2/table.hpp>
#include <nt2/include/constants/zero.hpp>
#include <nt2/include/constants/mone.hpp>
#include <nt2/include/functions/eps.hpp>
#include <nt2/include/functions/rec.hpp>
#include <nt2/include/functions/max.hpp>
#include <nt2/include/functions/min.hpp>
//TO DO !!
// #include <nt2/include/functions/diag.hpp>
// #include <nt2/include/functions/expand.hpp>
// #include <nt2/include/functions/fliplr.hpp>
// #include <nt2/include/functions/globalsum.hpp>
// #include <nt2/include/functions/trans.hpp>
// #include <nt2/include/functions/range.hpp>
// #include <nt2/include/functions/first_index.hpp>
// #include <nt2/include/functions/last_height_index.hpp>
#include <iostream>


//==============================================================================
// svd actual functor forward declaration
//==============================================================================
namespace nt2
{
  template<class A> struct svd_return;
  struct AllowDestroy {}; 
} 

namespace nt2 { namespace ext
{
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::svd_, tag::cpu_, 
                              (A)(S0)(B),
                              ((expr_< table_<unspecified_<A>,S0>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                              (scalar_<integer_<B> > )
                              )
  {
    typedef nt2::svd_return<A> result_type; 
    BOOST_DISPATCH_FORCE_INLINE result_type operator()(const A& a, const char & jobz) const
    {
      return nt2::svd_return<A>(a, jobz);
    }
  };  
} }

namespace nt2 
{
  //============================================================================
  // svd actual functor : precompute
  //============================================================================
  template<class A > struct svd_return
  {
    typedef typename A::value_type                   type_t;
    typedef typename A::index_type                  index_t; 
    typedef typename meta::as_real<type_t>::type    btype_t; 
    typedef nt2::table<type_t, nt2::matlab_index_>   ctab_t;
    typedef nt2::table<btype_t, nt2::matlab_index_> cbtab_t;
    typedef nt2::table<type_t, index_t>               tab_t;
    typedef nt2::table<btype_t, index_t>             btab_t;
    
    template < class XPR > svd_return(const XPR& a_, char jobz_ = 'A'):
      jobz(jobz_),
      a(a_),
      ma(a), 
      m(size(a, 1)),
      n(size(a, 2)),
      lda(leading_size(a)), 
      wrk(inw)
    {
      allocate(); 
      init(); 
    }
    template < class XPR>
    svd_return(const XPR & a_,
          nt2::details::workspace < type_t > & w_,
          char jobz_ = 'A'):
      jobz(jobz_),
      a(a_),
      ma(a),
      m(size(a, 1)),
      n(size(a, 2)),
      lda(leading_size(a)), 
      wrk(w_)        
    {
      allocate(); 
      init(); 
    }
      
//     template <            > svd_return( A &a,
//                                  nt2::details::workspace < type_t > & w_, 
//                                  details::allowDestroy destroy, char jobz_ = 'A'):
//       jobz(jobz_),
//       ma(a),
//       m(size(a, 1)),
//       n(size(a, 2)),
//       lda(nt2::details::padding(boost::proto::value(a))),
//       wrk(w_)        
//     {
//       allocate(); 
//       init(); 
//     }
//     template <            > svd_return( A &a,
//                                  details::allowDestroy destroy,
//                                  char jobz_ = 'A'):
//       jobz(jobz_),
//       ma(a),
//       m(size(a, 1)),
//       n(size(a, 2)),
//       lda(nt2::details::padding(boost::proto::value(a))),
//       wrk(inw)
//     {
//       allocate(); 
//       init(); 
//     }
    
    ~svd_return(){}
    // /////////////////////////////////////////////////////////////////////////////
    // accessors
    // /////////////////////////////////////////////////////////////////////////////
    tab_t       get_u ()      const
    {
      BOOST_ASSERT_MSG(nt2::details::lower(jobz) != 'n', "please call svd wit jobz = 'A', 'S' or 'O'"); 
      return u;
    }
    tab_t       get_vt()      const {
      BOOST_ASSERT_MSG(nt2::details::lower(jobz) != 'n', "please call svd wit jobz = 'A', 'S' or 'O'");     
      return vt;
    }
    btab_t      get_singular()const { return w; }
    btab_t      get_w()       const { return w; /*nt2::expand(nt2::diag(w), ucol, height(vt));*/}
    long int    get_info()    const { return info; }
    
    // /////////////////////////////////////////////////////////////////////////////
    // properties
    // /////////////////////////////////////////////////////////////////////////////
    btype_t     cond()       const { return  w(0)/w(nt2::min(m, n)-1); }
    size_t      rank(btype_t epsi = -1) const
    {
      return 0; /*nt2::globalsum(w > (epsi < 0 ? nt2::max(n, m)*nt2::eps(w(1)): epsi));*/
    }
    
    
    // /////////////////////////////////////////////////////////////////////////////
    // resolutions
    // /////////////////////////////////////////////////////////////////////////////
//     template < class XPR >tab_t solve(const XPR & b, btype_t epsi = Mone<btype_t>() )const{
//       epsi =  epsi < 0 ? nt2::eps(w(1)) : epsi; 
//       tab_t w1 = if_else( (w > epsi), nt2::rec(w), Zero<btype_t>());
//       return (nt2::trans(vt)*(nt2::diag(w1)*nt2::trans(u)))*b; 
//       //        return prodtMM(vt, prodMtM(diag(w1), u))*b;
//       }
      
//       template < class XPR > void solveip(XPR & b, btype_t epsi = -1 )const{
//         epsi =  epsi < 0 ? nt2::eps(w(1)) : epsi; 
//         ctab_t w1 = get_w();
//         w1 = if_else( (w1 > epsi), nt2::rec(w1), Zero<btype_t>());
//         b =  (nt2::trans(vt)*(nt2::diag(w1)*nt2::trans(u)))*b; 
//         //        b = prodtMM(vt, prodMtM(w1, u))*b;
//       }
      
//       tab_t null(btype_t epsi = -1 )const
//       {
//         epsi =  epsi < 0 ? nt2::eps(w(1)) : epsi;
//         // TODO use a reverse iterator on w
//         int j = length(w); 
//         for(; (j > 0) && (w(j)<= epsi); j--);
//         j++;
//         return nt2::fliplr(nt2::trans(vt(_(j, End()), _)));
//       }
      
//       tab_t orth(btype_t epsi =  -1)const
//       {
//         return u(_, nt2::_(1, rank(epsi))); 
//       }
      
//       tab_t zerosolve()const
//       {
//         return nt2::trans(vt(vt.last_height_index(), _));
//       }
      
//       tab_t pinv(btype_t epsi = -1 )const
//       {
//         epsi = epsi < 0 ? nt2::eps(w(1)) : epsi; 
//         tab_t w1 = nt2::trans(get_w());
//         w1 = if_else( (w1 > length(a)*epsi), rec(w1), Zero<btype_t>());
//         return (nt2::trans(vt)*(nt2::diag(w1)*nt2::trans(u)));
//         //        return prodtMM(vt, prodMtM(w1, u)); 
//       }
      
    private :
      inline void allocate()
      {
        switch (jobz)
          {
          case 'A':
            ldu = m;  ucol = m; 
            ldvt = n; vtcol = n; 
            break;
          case 'S':
            ldu = m;  ucol = nt2::min(n, m); 
            ldvt = n; vtcol = n; 
            break;
          case 'O':
            ldu = m;  ucol = nt2::min(n, m); 
            ldvt = n; vtcol = nt2::min(n, m);             
            break;
          case 'N':
            ldu = 1;  ucol = 1; 
            ldvt = 1; vtcol = 1;             
            break;
          default :
            break;
          }
        u.resize(of_size(ldu, ucol));
        ldu = leading_size(u); 
        vt.resize(of_size(ldvt, vtcol));
        ldvt = leading_size(vt); 
        w.resize(of_size(nt2::max(n, m), 1)); 
      }
      inline void init()
      {
        
        // gesdd se comporte moins bien que gesvd sur la matrice de hilbert 100x100
        // et plante pour jobz =  'N' (???JTL)
        //         gesdd(&jobz, &m, &n, ma.begin(), &lda, w.begin(),
        //               u.begin(), &ldu, vt.begin(), &ldvt,
        //               &info, wrk);
        
        // gesvd  ne donne pas exactement les même résultat que matlab 7.2 (???JTL) sur hilb(100)
        // légèrement moins bons, la recnstruction est à 8e-16 au lieu de 3e-16
        nt2::details::gesvd(&jobz, &jobz, &m, &n, ma.raw(), &lda,
                             w.raw(), u.raw(), &ldu,
                             vt.raw(), &ldvt, &info, wrk);
        
         for(size_t i = nt2::min(n, m)+1; i<= nt2::numel(w); i++) w(i) = Zero<btype_t>();
        //        mc_t::LapackTest(__FILE__, __LINE__, "gesvd", ma, info); 
      }

    const char                         jobz;
    A                                     a;
    A&                                   ma; 
    const long int                     m, n;
    const long int                      lda; 
    long int         ldu, ucol, ldvt, vtcol; 
    ctab_t                                u;
    ctab_t                               vt;
    cbtab_t                               w; 
    long int                           info; 
    nt2::details::workspace < type_t >  inw; 
    nt2::details::workspace < type_t > &wrk;
  };
} 
#endif
