/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_DETAILS_QR_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_DETAILS_QR_HPP_INCLUDED
// Matlab syntaxes
// Note that the syntaxes can be obtained using qr and pqr
// These 2 are qr ie no_p
//    [q,r] = qr(a), where a is m-by-n, produces an m-by-n upper triangular
//     matrix r and an m-by-m unitary matrix q so that a = q*r.
 
//     [q,r] = qr(a,0) produces the "economy size" decomposition.
//     if m>n, only the first n columns of q and the first n rows of r are
//     computed. if m<=n, this is the same as [q,r] = qr(a).

// and also these
//     x = qr(a) and x = qr(a,0) return the output of lapack's *geqrf routine.
//     (values() in the class)
//     triu(x) is the upper triangular factor r.
//////////// TODO
// to get the econmy size we much expand the results to the liked dimensions
////////////
// The others use pqr

//     if a is full:
 
//     [q,r,e] = qr(a) produces unitary q, upper triangular r and a
//     permutation matrix e so that a*e = q*r. the column permutation e is
//     chosen so that abs(diag(r)) is decreasing.
 
//     [q,r,e] = qr(a,'vector') returns the permutation information as a
//     vector instead of a matrix.  that is, e is a row vector such that 
//     a(:,e) = q*r. similarly, [q,r,e] = qr(a,'matrix') returns a permutation 
//     matrix e. this is the default behavior.
 
//     [q,r,e] = qr(a,0) produces an "economy size" decomposition in which e
//     is a permutation vector, so that a(:,e) = q*r.
 
 
#include <nt2/toolbox/linalg/details/utility/workspace.hpp>
#include <nt2/toolbox/linalg/details/utility/options.hpp>
#include <nt2/include/functions/qr.hpp>
#include <nt2/include/functions/of_size.hpp>
#include <nt2/include/functions/min.hpp>
#include <nt2/include/functions/max.hpp>
#include <nt2/include/functions/zeros.hpp>
#include <nt2/include/functions/eye.hpp>
#include <nt2/include/functions/triu.hpp>
#include <nt2/include/functions/expand.hpp>
#include <nt2/include/functions/prod.hpp>
#include <nt2/include/functions/height.hpp>
#include <nt2/include/functions/width.hpp>
//#include <nt2/include/functions/diag_of.hpp>
#include <nt2/include/constants/one.hpp>
#include <nt2/include/constants/eps.hpp>
#include <nt2/include/functions/isempty.hpp>
#include <nt2/toolbox/linalg/details/lapack/geqp3.hpp>
#include <nt2/toolbox/linalg/details/lapack/geqrf.hpp>
#include <nt2/toolbox/linalg/details/lapack/gqr.hpp>
#include <nt2/toolbox/linalg/details/lapack/mqr.hpp>
#include <nt2/toolbox/linalg/details/lapack/trtrs.hpp>
#include <nt2/table.hpp>
//#include <iostream>

namespace nt2 {
  struct no_p {};
  
  namespace details
  {
    template<class T> struct qr_result
    {
      typedef typename meta::strip<T>::type                   source_t;
      typedef typename source_t::value_type                     type_t;
      typedef typename meta::as_integer<type_t, signed>::type  itype_t;
      typedef typename source_t::index_type                    index_t;
      typedef typename meta::as_real<type_t>::type              base_t;
      typedef T                                                 data_t;
      typedef nt2::table<type_t,nt2::matlab_index_>              tab_t;
      typedef nt2::table<base_t,nt2::matlab_index_>             btab_t;
      typedef nt2::table<itype_t,nt2::matlab_index_>            itab_t;
      typedef nt2::details::workspace<type_t>              workspace_t;
      typedef nt2::table<nt2_la_int,nt2::matlab_index_>         ibuf_t;
      typedef nt2::table<type_t,index_t>                   result_type;
      
      template<class Input>
      qr_result ( Input& xpr)
        : a_(xpr)
        , aa_(xpr)
        , m_(nt2::height(xpr))
        , n_(nt2::width(xpr))
        , k_(nt2::min(m_, n_))
        , lda_(a_.leading_size())
        , q_(of_size(0, 1))
        , jpvt_(nt2::of_size(n_,1)) 
        , tau_(nt2::of_size(k_,1)) 
        , info_(0)
        , p_(of_size(0, 1))
        , nop_(false)
      {
        nt2::details::geqp3(&m_, &n_, aa_.raw(), &lda_,
                            jpvt_.raw(), tau_.raw(), &info_);
      }
      
      template<class Input>
      qr_result ( Input& xpr, const no_p&)
        : a_(xpr)
        , aa_(xpr)
        , m_(nt2::height(xpr))
        , n_(nt2::width(xpr))
        , k_(nt2::min(m_, n_))
        , lda_(a_.leading_size())
        , q_(of_size(0, 1))
        , jpvt_(nt2::of_size(0,1)) 
        , tau_(nt2::of_size(k_,1)) 
        , info_(0)
        , p_(of_size(0, 1))
        , nop_(true)
      {
        nt2::details::geqrf(&m_, &n_, aa_.raw(), &lda_, tau_.raw(), &info_, wrk_);
      }
      
      qr_result& operator=(qr_result const& src)
      {
        a_      = src.a_;
        aa_     = src.aa_;
        m_      = src.m_;
        n_      = src.n_;
        k_      = src.k_;
        lda_    = src.lda_;
        jpvt_   = src.jpvt_; 
        tau_    = src.tau_;  
        info_   = src.info_;
        p_ = src.p_;
        q_ = src.q_;
        nop_ =  src.nop_; 
        return *this;
      }
      result_type values() const { return aa_; }
      
      result_type q ()
      {
        if(is_empty(q_))
          {
            nt2_la_int nn = nop_ ? k_ : m_; 
            q_ = nt2::expand(aa_, nn, nn);
            nt2::details::gqr(&m_, &nn, &k_, q_.raw(), &lda_, tau_.raw(), &info_);
          }
        return q_;
      }
      result_type r()const
      {
        return triu(aa_);
      }
      result_type p()
      {
        if(is_empty(p_))
          {
            if (nop_)
              {
                p_ = nt2::eye(n_, n_, meta::as_<type_t>());
                return p_; 
              }
            p_ = nt2::zeros(nt2::numel(jpvt_), nt2::meta::as_<type_t>()); 
            for(unsigned int i=1; i <= nt2::size(p_, 1) ; ++i){
              p_(jpvt_(i), i) = One<type_t>(); 
            }
          }
        return p_; 
      }
      
      itab_t  jp()
      {
        itab_t jip(of_size(1, numel(jpvt_))); 
        if (nop_)
          {
            return ones(1, n_, nt2::meta::as_<itype_t>());
          }
        for(size_t i=1; i <= numel(jpvt_); ++i) jip(i) = jpvt_(i); 
        return jip; 
      }
      
      nt2_la_int status() const { return info_; }
      
      size_t     rank(base_t epsi = nt2::Eps<base_t>())const
      {
        btab_t m = nt2::max(nt2::abs(diag_of(aa_))); 
        base_t thresh = nt2::max(n_, m_)*epsi*m(1);
        size_t r = 0; 
        for(int i=1; i <= min(n_, m_); ++i)
          {
            if (nt2::abs(aa_(i, i)) > thresh) ++r; 
          }
        return r;
        //nt2::sum(nt2::abs(nt2::diag_of(aa_())) > nt2::max(n, m)*epsi*nt2::max(nt2::abs(aa_)));
      }
      base_t absdet()const{
        BOOST_ASSERT_MSG(m_ == n_, "non square matrix in determinant computation"); 
        return nt2::prod(nt2::abs(diag_of(aa_)))(1);
      }
      
      //==========================================================================
      // Solver interface
      //==========================================================================
      template<class XPR> result_type solve(const XPR & b, base_t epsi = nt2::Eps<base_t>(),
                                            bool transpose = false)const
      {
        result_type bb = b;
        inplace_solve(bb, epsi, transpose);
        return bb;
      }
      
      template < class XPR > void inplace_solve(XPR & b, base_t epsi = nt2::Eps<base_t>(),
                                                bool transpose = false) const
      {
        char side = 'L';
        char tr = (transpose) ? 'N' : !is_real(type_t(1))? 'C':'T';
        nt2_la_int M = nt2::size(b, 1), N = nt2::size(b, 2); 
        nt2_la_int ldb = b.leading_size();
        nt2_la_int info; 
        nt2::details::mqr(&side, &tr, &M, &N, &k_, aa_.raw(), &lda_, tau_.raw(), b.raw(), &ldb, &info, wrk_);
        nt2_la_int nrhs = size(b, 2); 
        char uplo =  'U', d = 'N';
        char tr1 = (transpose) ?  (!is_real(type_t(1))? 'C':'T') : 'N';
        nt2_la_int rk = rank(epsi); 
        nt2::details::trtrs(&uplo, &tr1, &d, &rk, &nrhs, aa_.raw(), &lda_, b.raw(), &ldb, &info);
        if (!nop_) b = permute(b);
      }
      
    private :
      
      inline tab_t permute(const tab_t& bb) const {
        tab_t res(nt2::of_size(nt2::numel(jpvt_), nt2::size(bb, 2)));
        const size_t m =  nt2::min(size(bb, 1), numel(jpvt_));
        for(int i=1; i <= m; ++i)
          {
            res(jpvt_(i), _) = bb(i, _); 
          }
        //      res(jpvt_(_(1, m)), _) = bb; //TODO
        return res; 
      }
      template < class S>
      static S diag_of(const S& a)
      {
        S d(of_size(nt2::min(width(a), height(a)), 1)); 
        for (int i = 1; i <= nt2::min(width(a), height(a)); ++i) d(i) = a(i, i);
        return d; 
      }
      
      data_t                 a_; 
      tab_t                 aa_;
      nt2_la_int     m_, n_, k_;
      nt2_la_int           lda_; 
      tab_t                  q_; 
      ibuf_t              jpvt_;
      tab_t                tau_; 
      nt2_la_int          info_;
      tab_t                  p_;
      bool                 nop_;
      mutable workspace_t  wrk_;

    };
    
  }
}
#endif
