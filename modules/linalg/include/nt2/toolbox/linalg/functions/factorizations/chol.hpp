/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_FACTORIZATIONS_CHOL_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_FACTORIZATIONS_CHOL_HPP_INCLUDED

#include <nt2/include/functions/chol.hpp>
#include <nt2/include/functions/numel.hpp>
#include <nt2/toolbox/linalg/details/utility/padding.hpp>
#include <nt2/toolbox/linalg/details/lapack/chol.hpp>
#include <nt2/toolbox/linalg/details/lapack/con.hpp>
#include <nt2/toolbox/linalg/details/lapack/lange.hpp>

#include <nt2/table.hpp>
//#include <iostream>
//#include <nt2/include/functions/triu.hpp>


//==============================================================================
// chol actual functor forward declaration
//==============================================================================
namespace nt2
{
  template<class A> struct chol_f;
  struct AllowDestroy {}; 
} 

namespace nt2 { namespace ext
{
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::chol_, tag::cpu_, 
                              (A)(S0),
                              ((expr_< table_<unspecified_<A>,S0>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                              )
  {
    typedef nt2::chol_f<A> result_type; 
    BOOST_DISPATCH_FORCE_INLINE result_type operator()(const A& a) const
    {
      return nt2::chol_f<A>(a);
    }
  };  
} }

namespace nt2 
{
  //============================================================================
  // chol actual functor : precompute
  //============================================================================
  template<class A > struct chol_f
  {
    typedef typename A::value_type                   type_t;
    typedef typename A::index_type                  index_t; 
    typedef typename meta::as_real<type_t>::type     base_t; 
    typedef nt2::table<type_t, nt2::matlab_index_>   ftab_t;
    typedef nt2::table<type_t, nt2::matlab_index_>  fbtab_t;
    typedef nt2::table<type_t, index_t>               tab_t;
    typedef nt2::table<type_t, index_t>              btab_t;
    
      template < class XPR > chol_f(const XPR & a_):
        uplo('U'),
        a(a_),
        ma(a),
        n(size(a, 2)),
        lda(boost::proto::value(a).leading_size()) 
      {
        //assert issquare a_
        init(); 
      }
      template < class T, class S > chol_f( table<T, S>&a, AllowDestroy destroy):
        uplo('U'),
        ma(a),
        n(size(a, 2)),
        lda(nt2::details::padding(boost::proto::value(a)))
      {
        init(); 
      }
      
      ~chol_f(){}
      // /////////////////////////////////////////////////////////////////////////////
      // accessors
      // /////////////////////////////////////////////////////////////////////////////
    tab_t       getu () const   { return ma; }//triu(ma); }
      long int    getinfo() const { return info;     }


      // /////////////////////////////////////////////////////////////////////////////
      // computations
      // /////////////////////////////////////////////////////////////////////////////
      size_t negativity(){ return (info < 0) ? 0 : info; }
      base_t rcond()
      {
        char norm = '1';
        base_t res;
        base_t anorm = nt2::details::lange(&norm,  &n,  &n, a.raw(), &lda); 
        nt2::details::pocon(&uplo, &n, a.raw(), &lda, &anorm, &res, &info); 
        return res;  
      }
     
      // /////////////////////////////////////////////////////////////////////////////
      // resolutions
      // /////////////////////////////////////////////////////////////////////////////
      template < class XPR > ftab_t solve(const XPR & b ){
        ftab_t bb = b; 
        long int nrhs = size(bb, 1);
        long int ldb  = leading_size(bb);
        nt2::details::potrs(&uplo, &n, &nrhs, ma.raw(), &lda, bb.raw(), &ldb, &info);
        return bb; 
      }
      
      template < class XPR > void solveip(XPR & b ){
        long int nrhs = size(b, 1);
        long int ldb  = leading_size(b);
        nt2::details::potrs(&uplo, &n, &nrhs, ma.raw(), &lda, b.raw(), &ldb, &info);
     }
      
      
    private :

      inline void init()
      {
        nt2::details::potrf(&uplo, &n, ma.raw(), &lda, &info);
        //        mc_t::LapackTest(__FILE__, __LINE__, "potrf", ma, info); 
      }
      const char     uplo;
      A                 a;
      A&               ma; //ma has to be a view
      const long int    n;
      const long int  lda; 
      long int       info;
    };
} 
#endif
