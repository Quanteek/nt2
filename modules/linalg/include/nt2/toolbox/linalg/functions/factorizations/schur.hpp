/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_FACTORIZATIONS_SCHUR_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_FACTORIZATIONS_SCHUR_HPP_INCLUDED

#include <nt2/toolbox/linalg/details/lapack/workspace.hpp>
#include <nt2/toolbox/linalg/details/utility/options.hpp>
#include <nt2/toolbox/linalg/details/utility/lower.hpp>
#include <nt2/toolbox/linalg/details/utility/tags.hpp>
#include <nt2/include/constants/eps.hpp>
#include <nt2/include/functions/schur.hpp>
#include <nt2/include/functions/of_size.hpp>
#include <nt2/include/functions/min.hpp>
#include <nt2/include/functions/isempty.hpp>
#include <nt2/include/functions/leading_size.hpp>
#include <nt2/toolbox/linalg/details/lapack/geesx.hpp>

#include <nt2/table.hpp>
//#include <iostream>
//#include <nt2/include/functions/expand.hpp>
//#include <nt2/include/functions/diag.hpp>
//#include <nt2/include/functions/prod.hpp>
//#include <nt2/include/functions/triu.hpp>
//#include <nt2/include/functions/range.hpp>


//==============================================================================
// schur call result forward declaration
//==============================================================================
namespace nt2
{
  template<class A> struct schur_return;

} 

namespace nt2 { namespace ext
{
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::schur_, tag::cpu_, 
                              (A)(S0),
                              ((expr_< table_<unspecified_<A>,S0>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                              )
  {
    typedef nt2::schur_return<A> result_type; 
    BOOST_DISPATCH_FORCE_INLINE result_type operator()(const A& a) const
    {
      return nt2::schur_return<A>(a);
    }
  };  
} }

namespace nt2 
{
  //============================================================================
  // schur actual functor : precompute
  //============================================================================
  template<class A,
           class CPLX = typename nt2::meta::is_complx<typename A::value_type >::type>
  struct schur_return
  {
    typedef typename A::value_type                     type_t;
    typedef typename A::index_type                    index_t; 
    typedef typename meta::as_real<type_t>::type       base_t; 
    typedef nt2::table<type_t, nt2::matlab_index_>     ftab_t;
    typedef nt2::table<type_t, nt2::matlab_index_>    fbtab_t;
    typedef nt2::table<long int, nt2::matlab_index_>  iftab_t;
    typedef nt2::table<type_t, index_t>                 tab_t;
    typedef nt2::table<long int, index_t>              itab_t;
    typedef nt2::table<type_t, index_t>                btab_t;
    typedef nt2::details::workspace<type_t>       workspace_t;
    typedef nt2::details::options                   options_t;
    typedef nt2::details::allowDestroy         allowDestroy_t;
    typedef std::complex<base_t>                      ctype_t; 
    typedef nt2::table<ctype_t, nt2::matlab_index_>   cftab_t;
    typedef nt2::table<ctype_t, nt2::matlab_index_>  cfbtab_t;
    typedef nt2::table<ctype_t, index_t>              cftab_t;
    typedef nt2::table<ctype_t, index_t>             cfbtab_t;

     
    static inline long int selectall(void* a) { return true; }
    typedef long int select_t(base_t*); 
     
      template < class XPR > schur_return(const XPR & a_,
                                   char jobvs_ = 'V',
                                   char sort_  = 'N',
                                   char sense_ = 'N'
                                   ):
        jobvs(nt2::details::lower(jobvs_)),
        sort(nt2::details::lower(sort_)), 
        sense(nt2::details::lower(sense_)), 
        a(a_),
        ma(a),
        n(nt2::size(a, 1)),
        lda(leading_size(a),
        wrk(inw)
      { 
        allocate(); 
        init();
      }
      template < class XPR > schur_return(const XPR & a_,
                                   workspace_t & w_, 
                                   char jobvs_ = 'V', 
                                   char sort_  = 'N',
                                   char sense_ = 'N'
                                   ):
        jobvs(nt2::details::lower(jobvs_)),
        sort( nt2::details::lower(sort_)), 
        sense(nt2::details::lower(sense_)), 
        ma(a),
        n(nt2::size(a, 1)),
        lda(leading_size(a),
        wrk(w_)
      {
        allocate(); 
        init(); 
      }
      template < class XPR > schur_return(const XPR & a_,
                                   allowDestroy_t destroy, 
                                   char jobvs_ = 'V',
                                   char sort_  = 'N', 
                                   char sense_ = 'N'
                                   ):
        jobvs(nt2::details::lower(jobvs_)), 
        sort( nt2::details::lower(sort_)), 
        sense(nt2::details::lower(sense_)), 
        ma(a),
        n(nt2::size(a, 1)),
        lda(leading_size(a),
        wrk(inw)
      {
        allocate(); 
        init(); 
      }
      template < class XPR > schur_return(const ne::expression<XPR> & a_,
                                   workspace_t & w_, 
                                   allowDestroy_t destroy, 
                                   char jobvs_ = 'V',
                                   char sort_  = 'N', 
                                   char sense_ = 'N'
                                   ):
        jobvs(nt2::details::lower(jobvs_)), 
        sort( nt2::details::lower(sort_)), 
        sense(nt2::details::lower(sense_)), 
        a(a_),
        ma(a),
        n(nt2::size(a, 1)),
        lda(leading_size(a),
        wrk(w_)
      {
        allocate(); 
        init(); 
      }
       
        
      ~schur_return(){}
      // /////////////////////////////////////////////////////////////////////////////
      // accessors
      // /////////////////////////////////////////////////////////////////////////////
      //      tab_t       getvl  ()        { return vl;     }
      //      tab_t       getvr  ()        { return vr;     }
      //      ctab_t      gete   ()        { return w;      }
      ctab_t       getw   () const     { return diag(w);}
      ctab_t       gett   ()           { return ma;     }
      ctab_t       getz   ()
      {
        //       mch_t::LapackOption(__FILE__, __LINE__, "jobs", jobvs == 'v', jobvs); 
        return vs;
      }
      long int    getinfo() const     { return info;   }
      base_t      getabnorm() const   { return abnorm; }
      base_t      getrconde()
      {
        //        mch_t::LapackOption(__FILE__, __LINE__, "sense", (sense ==  'e' || sense ==  'b'), sense); 
        return rconde; 
      }
      base_t      getrcondv()
      {
        //        mch_t::LapackOption(__FILE__, __LINE__, "sense", (sense ==  'e' || sense ==  'b'), sense); 
        return rcondv; 
      }
     
    private :
      inline void allocate() 
      {
        jobvs = (sense == 'e' || sense == 'b') ? 'v':jobvs;
        sort = (sense == 'n')? sort : 's'; 
        ldvs = (jobvs == 'v') ? n : 1;
        w.resize(nt2::of_size(n, 1));
        vs.resize(of_size(ldvs, ldvs));
        ldvs = leading_size(vs); 
      }
      
      inline void init()
      {
        nt2::details::geesx(&jobvs, &sort, &selectall , &sense, &n, 
             ma.raw(), &lda, &sdim, w.raw(),
             vs.raw(), &ldvs,
             &rconde, &rcondv, 
             &info, wrk);
        //       mch_t::LapackTest(__FILE__, __LINE__, "geesx", ma, info); 
      }
      char                 jobvs, sort;
      const char                 sense;
      ftab_t                         a;
      ftab_t                        ma; 
      const long int                 n;
      const long int               lda; 
      long int                    ldvs;
      long int                    sdim; 
      base_t                    abnorm; 
      base_t            rcondv, rconde;
      cftab_t                       vs; 
      cfta_t                         w; 
      long int                    info; 
      workspace_t                  inw; 
      workspace_t                 &wrk;
    };

    // //////////////////////////////////////////////////////////////////////////////////////
    // real case
    // //////////////////////////////////////////////////////////////////////////////////////
  template<class A > 
            struct schur_return < A, boost::mpl_::false_ > 
  {
    typedef typename A::value_type                     type_t;
    typedef typename A::index_type                    index_t; 
    typedef typename meta::as_real<type_t>::type       base_t; 
    typedef nt2::table<type_t, nt2::matlab_index_>     ftab_t;
    typedef nt2::table<type_t, nt2::matlab_index_>    fbtab_t;
    typedef nt2::table<long int, nt2::matlab_index_>  iftab_t;
    typedef nt2::table<type_t, index_t>                 tab_t;
    typedef nt2::table<long int, index_t>              itab_t;
    typedef nt2::table<type_t, index_t>                btab_t;
    typedef nt2::details::workspace<type_t>       workspace_t;
    typedef nt2::details::options                   options_t;
    typedef nt2::details::allowDestroy         allowDestroy_t;
    typedef std::complex<base_t>                      ctype_t; 
    typedef nt2::table<ctype_t, nt2::matlab_index_>   cftab_t;
    typedef nt2::table<ctype_t, nt2::matlab_index_>  cfbtab_t;
    typedef nt2::table<ctype_t, index_t>              cftab_t;
    typedef nt2::table<ctype_t, index_t>             cfbtab_t;

    static inline long int selectall2(base_t* a, base_t* b) { return true; }
     
    typedef long int select2_t(base_t*, base_t*); 
     
      template < class XPR > schur_return(const XPR & a_,
                                   char jobvs_ = 'V',
                                   char sort_  = 'N',
                                   char sense_ = 'N'
                                   ):
        jobvs(nt2::details::lower(jobvs_)),
        sort(nt2::details::lower(sort_)), 
        sense(nt2::details::lower(sense_)), 
        a(a_),
        ma(a),
        n(nt2::size(a, 1)),
        lda(leading_size(a),
        wrk(inw)
      { 
        allocate(); 
        init();
      }
      template < class XPR > schur_return(const XPR & a_,
                                   workspace_t & w_, 
                                   char jobvs_ = 'V', 
                                   char sort_  = 'N',
                                   char sense_ = 'N'
                                   ):
        jobvs(nt2::details::lower(jobvs_)),
        sort( nt2::details::lower(sort_)), 
        sense(nt2::details::lower(sense_)), 
        ma(a),
        n(nt2::size(a, 1)),
        lda(leading_size(a),
        wrk(w_)
      {
        allocate(); 
        init(); 
      }
      template < class XPR > schur_return(const XPR & a_,
                                   allowDestroy_t destroy, 
                                   char jobvs_ = 'V',
                                   char sort_  = 'N', 
                                   char sense_ = 'N'
                                   ):
        jobvs(nt2::details::lower(jobvs_)), 
        sort( nt2::details::lower(sort_)), 
        sense(nt2::details::lower(sense_)), 
        ma(a),
        n(nt2::size(a, 1)),
        lda(leading_size(a),
        wrk(inw)
      {
        allocate(); 
        init(); 
      }
      template < class XPR > schur_return(const ne::expression<XPR> & a_,
                                   workspace_t & w_, 
                                   allowDestroy_t destroy, 
                                   char jobvs_ = 'V',
                                   char sort_  = 'N', 
                                   char sense_ = 'N'
                                   ):
        jobvs(nt2::details::lower(jobvs_)), 
        sort( nt2::details::lower(sort_)), 
        sense(nt2::details::lower(sense_)), 
        a(a_),
        ma(a),
        n(nt2::size(a, 1)),
        lda(leading_size(a),
        wrk(w_)
      {
        allocate(); 
        init(); 
      }
       
        
      ~schur_return(){}
      // /////////////////////////////////////////////////////////////////////////////
      // accessors
      // /////////////////////////////////////////////////////////////////////////////
      ctab_t        getw   () const     { return diag(gete());}
      ctab_t        gete   () const     { return cmplx(wr, wi); }
      tab_t         geter  ()           { return wr;     }
      tab_t         getei  ()           { return wi;     }   
      tab_t         gett   ()
      {
        return ma;
      }
      tab_t       getz   ()
      {
        //      mch_t::LapackOption(__FILE__, __LINE__, "jobs", jobvs == 'v', jobvs); 
        return vs;
      }
      long int    getinfo() const     { return info;   }
      base_t      getabnorm() const   { return abnorm; }
      base_t      getrconde()
      {
        //        mch_t::LapackOption(__FILE__, __LINE__, "sense", (sense ==  'e' || sense ==  'b'), sense); 
        return rconde; 
      }
      base_t      getrcondv()
      {
        //        mch_t::LapackOption(__FILE__, __LINE__, "sense", (sense ==  'e' || sense ==  'b'), sense); 
        return rcondv; 
      }
      
    private :
      inline void allocate() 
      {
        jobvs = (sense == 'e' || sense == 'b') ? 'v':jobvs;
        sort = (sense == 'n')? sort : 's'; 
        ldvs = (jobvs == 'v') ? n : 1;
        wr.resize(nt2::of_size(n, 1));
        wi.resize(nt2::of_size(n, 1)); 
        vs.resize(nt2::of_size(ldvs, ldvs));
        ldvs = leading_size(vs); 
      }
      
      inline void init()
      {
        nt2::details::geesx(&jobvs, &sort, &selectall2, &sense, &n, 
              ma.raw(), &lda, &sdim, wr.raw(), wi.raw(), 
              vs.raw(), &ldvs,
              &rconde, &rcondv, 
              &info, wrk);
        //     mch_t::LapackTest(__FILE__, __LINE__, "geesx", ma, info); 
      }
     char                  jobvs, sort;
      const char                 sense;
      ftab_t                         a;
      ftab_t                        ma; 
      const long int                 n;
      const long int               lda; 
      long int                    ldvs;
      long int                    sdim; 
      base_t                    abnorm; 
      base_t            rcondv, rconde;
      cftab_t                       vs; 
      ftab_t                    wr, wi; 
      long int                    info; 
      workspace_t                  inw; 
      workspace_t                 &wrk;
  };
} 
#endif
