/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_FACTORIZATIONS_PLU_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_FACTORIZATIONS_PLU_HPP_INCLUDED

#include <nt2/toolbox/linalg/details/lapack/workspace.hpp>
#include <nt2/toolbox/linalg/details/utility/options.hpp>
#include <nt2/include/constants/eps.hpp>
#include <nt2/include/functions/plu.hpp>
#include <nt2/include/functions/of_size.hpp>
#include <nt2/include/functions/min.hpp>
#include <nt2/include/functions/isempty.hpp>
#include <nt2/include/functions/leading_size.hpp>
#include <nt2/toolbox/linalg/details/lapack/lange.hpp>
#include <nt2/toolbox/linalg/details/lapack/gecon.hpp>
#include <nt2/toolbox/linalg/details/lapack/getri.hpp>
#include <nt2/toolbox/linalg/details/lapack/gesvx.hpp>
#include <nt2/toolbox/linalg/details/lapack/getrf.hpp>
#include <nt2/include/functions/triu.hpp>
#include <nt2/include/functions/tril.hpp>
#include <nt2/include/functions/tri1l.hpp>
#include <nt2/include/functions/eye.hpp>

#include <nt2/table.hpp>
//#include <iostream>
//#include <nt2/include/functions/expand.hpp>
//#include <nt2/include/functions/diag.hpp>
//#include <nt2/include/functions/prod.hpp>
//#include <nt2/include/functions/range.hpp>

// rref computation
// The U matrix in an LU decomposition is a ref. Dividing each row by the leading entry will convert a ref to a rref. 

// plu rank is in status

// To get the sign of the determinant of the permutation matrix, you start with + and you go 
// through the IPIV array for I=1 to N and you change sign each time IPIV(I) is not I. (You might 
// want stop at N-1 and not N since IPIV(N)=N.) Then you product the diagonal term of U and 
// multiply with the sign of the determinant of P and you obtain the determinant of A from the 
// A=PLU factorization. 

// There is determinant computation in LINPACK, see for example the routine DGEDI. 
// http://www.netlib.org/linpack/dgedi.f 
// Computing a determinant is likely to overflow, the LINPACK's routine is specially careful 
// about that.
// c        det     double precision(2)
// c                determinant of original matrix if requested.
// c                otherwise not referenced.
// c                determinant = det(1) * 10.0**det(2)
// c                with  1.0 .le. dabs(det(1)) .lt. 10.0
// c                or  det(1) .eq. 0.0 .
// c     compute determinant
// c
// TODO perhaps a version using ldexp frexp 2^n rather that 10^n
//       if (job/10 .eq. 0) go to 70
//          det(1) = 1.0d0
//          det(2) = 0.0d0
//          ten = 10.0d0
//          do 50 i = 1, n
//             if (ipvt(i) .ne. i) det(1) = -det(1)
//             det(1) = a(i,i)*det(1)
// c        ...exit
//             if (det(1) .eq. 0.0d0) go to 60
//    10       if (dabs(det(1)) .ge. 1.0d0) go to 20
//                det(1) = ten*det(1)
//                det(2) = det(2) - 1.0d0
//             go to 10
//    20       continue
//    30       if (dabs(det(1)) .lt. ten) go to 40
//                det(1) = det(1)/ten
//                det(2) = det(2) + 1.0d0
//             go to 30
//    40       continue
//    50    continue
//    60    continue
//    70 continue


//==============================================================================
// plu call result forward declaration
//==============================================================================
namespace nt2
{
  template<class A> struct plu_return;

} 

namespace nt2 { namespace ext
{
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::plu_, tag::cpu_, 
                              (A)(S0),
                              ((expr_< table_<unspecified_<A>,S0>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                              )
  {
    typedef nt2::plu_return<A> result_type; 
    BOOST_DISPATCH_FORCE_INLINE result_type operator()(const A& a) const
    {
      return nt2::plu_return<A>(a);
    }
  };  
} }

namespace nt2 
{
  //============================================================================
  // plu actual functor : precompute
  //============================================================================
  template<class A > struct plu_return
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
    
    template < class XPR > plu_return(const XPR & a_):
      a(a_), 
      lu(a_),
      m(size(a, 1)),
      n(size(a, 2)),
      lda(leading_size(a)), 
      ipiv(nt2::of_size(nt2::min(n, m), 1)), 
      r(of_size(n, 1)),
      c(of_size(n, 1)), 
      rc(-1), 
      w(inw), 
      p_computed(false),
      pt_computed(false)
    {
      init();
    }
    template < class XPR > plu_return(const XPR & a_,
                                      workspace_t &w_):
      a(a_), 
      lu(a_),
      m(size(a, 1)),
      n(size(a, 2)),
      lda(leading_size(a)), 
      ipiv(nt2::of_size(nt2::min(n, m), 1)), 
      r(of_size(n, 1)),
      c(of_size(n, 1)), 
      rc(-1), 
      w(w_), 
      p_computed(false),
      pt_computed(false)
    {
      init();
    }
    
    ~plu_return(){}
    // /////////////////////////////////////////////////////////////////////////////
    // accessors
    // /////////////////////////////////////////////////////////////////////////////
    tab_t     geta()  {return a; }
    tab_t     getu()  {return triu(lu/*(_(1, std::min(n, m)),_) */);  }
    tab_t     getl()  {return  tri1l(lu/*tri1l(lu(_,_(1, std::min(n, m)))*/ ); }
    itab_t    getip() {return ipiv; }
    itab_t    getp()  { //TODO optimize the call
      if (pt_computed) return trans(p); 
      if (!p_computed)
        {
          p = nt2::eye(m, m, meta::as_<base_t>());
          for(size_t i=1; i <= numel(ipiv); ++i){
            tab_t c = p(i, _);
            p(i,_) = p(ipiv(i),_);
            p(ipiv(i),_) = c; 
          }
          p_computed = true; 
          return p; 
        }
      else
        {
          return p; 
        }
    } 
    tab_t     getpt(){//TODO optimize the call
      if (p_computed) return trans(p); 
      if (!pt_computed)
        {
          p = nt2::eye(m, m, meta::as_<base_t>());
          for(size_t i=1; i <= numel(ipiv); ++i){
            tab_t c = p(_,i);
            p(_,i) = p(_,ipiv(i));
            p(_,ipiv(i)) = c; 
          }
          pt_computed = true; 
          return p; 
        }
      else
        {
          return p;
        }
    } 
    
    long int getinfo(){return info; }
    
    // /////////////////////////////////////////////////////////////////////////////
    // computations
    // /////////////////////////////////////////////////////////////////////////////
    base_t rcond(char c = '1')
    {
      if (c !=  '1' || rc == -1)
        {
          char norm = c;
          base_t anorm = nt2::details::lange(&norm,  &n,  &n, lu.raw(), &lda, w); 
          nt2::details::gecon(&norm, &n,  lu.raw(), &lda, &anorm, &rc, &info, w);
        }
      return rc;  
    }
    
    
    size_t rank(base_t epsi = nt2::Eps<base_t>())
    {
      return 0; //globalSum( abs(diag(getu())) > std::max(n, m)*epsi*globalMax(abs(diag(getu()))) );
    }
    
    //       base_t absdet(){
    //         //        mc_t::SquareTest(__FILE__, __LINE__, a);
    //         return globalProd(abs(diag(getu())));
    //       }
    
    //       type_t det(){
    //         //       mc_t::SquareTest(__FILE__, __LINE__, a);
    //         return globalProd(diag(getu()))*((globalSum(ipiv != ne::convert<long>(colvect((colon(1, numel(ipiv))))))%2 == 1) ? -1 : 1); 
    //       }
    
    tab_t inv(bool noWarn =  false)
    {
      tab_t i = lu;
      base_t rc; 
      if (!noWarn && (rc = rcond()) < nt2::Eps<base_t>())
        {
          std::cerr << " Warning : Na::Matrix is close to singular or badly scaled." << std::endl; 
          std::cerr << "           Results may be inaccurate. RCOND = " << rc << "." << std::endl;
        }
      nt2::details::getri(&n, i.raw(), &lda, ipiv.raw(), &info, w);
      return i; 
    }  
    
    // /////////////////////////////////////////////////////////////////////////////
    // resolutions
    // /////////////////////////////////////////////////////////////////////////////
    template < class XPR >
    tab_t solve(const XPR & b, bool noWarn =  false)
    {
      tab_t x(size(b)), bb(b); 
      long int nrhs =  size(b, 2);
      long int ldb  =  leading_size(b);
      long int ldx  =  leading_size(x); 
      long int lda  =  leading_size(a);
      long int ldlu =  leading_size(lu);  
      if (isempty(berr))
        {
          berr.resize(of_size(nrhs, 1));
          ferr.resize(of_size(nrhs, 1));
        }
      nt2::details::gesvx((char*)&options_t::opt('F'),
                          (char*)&options_t::opt('N'),
                          &n, &nrhs,
                          a.raw(), &lda,
                          lu.raw(), &ldlu,
                          ipiv.raw(),
                          (char*)&options_t::opt('N'), 
                          r.raw(), c.raw(),
                          bb.raw(), &ldb,
                          x.raw(), &ldx,
                          &rc,
                          ferr.raw(),
                          berr.raw(),
                          &info, w);
      if(!noWarn && (info > n || rc < nt2::Eps<base_t>())){
        std::cerr << " Warning : Na::Matrix is close to singular or badly scaled." << std::endl; 
        std::cerr << "           Results may be inaccurate. RCOND = " << rc << "." << std::endl;
      }
      //        mc_t::LapackTest(__FILE__, __LINE__, "gesvx", a, info); 
      return x;
    }

    base_t rcond()
    {
      char norm = '1';
      base_t res;
      base_t anorm = nt2::details::lange(&norm,  &n,  &n, a.raw(), &lda); 
      nt2::details::gecon(&norm, &n, lu.raw(), &lda, &anorm, &res, &info); 
      return res;  
    }
    
  private :
    inline void init()
    {
      nt2::details::getrf(&m, &n, lu.raw(), &lda, ipiv.raw(), &info, w);
      //        mc_t::LapackTest(__FILE__, __LINE__, "getrf", lu, info); 
    }
    tab_t                 a,lu;
    const long int         m,n;
    const long int         lda; 
    itab_t                ipiv;
    tab_t                    u;
    tab_t                   vt;
    btab_t       r,c,ferr,berr; 
    long int              info; 
    btab_t                   p; 
    base_t                  rc;
    workspace_t            inw; 
    workspace_t             &w;
    bool            p_computed;
    bool           pt_computed; 
  };
} 
#endif
