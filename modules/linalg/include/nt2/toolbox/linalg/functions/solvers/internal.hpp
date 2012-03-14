#ifndef __PROJECT__FILE__NT2_ALGEBRA_SOLVERS_INTERNAL_HPP__INCLUDED
#define __PROJECT__FILE__NT2_ALGEBRA_SOLVERS_INTERNAL_HPP__INCLUDED

////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2007 for LASMEA UMR 6602 du CNRS.                  
//  All rights reserved.                                             
//                                                                   
//  License informations are available in the LICENSE file.          
//  Additionnal informations are available on http://nt2.sourceforge.net 
////////////////////////////////////////////////////////////////////////////////
#include <nt2/toolbox/linalg/details/lapack/gesv.hpp>
#include <nt2/toolbox/linalg/details/lapack/gels.hpp>
#include <nt2/toolbox/linalg/details/lapack/gelsd.hpp>  
#include <nt2/toolbox/linalg/details/lapack/posv.hpp>
#include <nt2/toolbox/linalg/details/lapack/trtrs.hpp>
#include <nt2/include/functions/height.hpp>
#include <nt2/include/functions/width.hpp>
#include <nt2/include/functions/leading_size.hpp>
#include <nt2/include/functions/issquare.hpp>
#include <nt2/include/functions/areofsameheight.hpp>
#include <nt2/include/functions/max.hpp>
#include <nt2/include/functions/min.hpp>
#include <nt2/sdk/meta/strip.hpp>
#include <nt2/table.hpp>
//#include <nt2/include/functions/expand.hpp>
//#include <nt2/include/functions/range.hpp>

namespace nt2
{
  namespace details
  {
    
    template < class XPR0,  class XPR1,  class XPR2 = XPR1> 
    class solvers
    {
    public:
    typedef long int                                 la_int; 
    typedef typename A::value_type                   type_t;
    typedef typename A::index_type                  index_t; 
    typedef typename meta::as_real<type_t>::type    btype_t; 
    typedef nt2::table<type_t,nt2::matlab_index_>    ftab_t;
    typedef nt2::table<btype_t,nt2::matlab_index_>  fbtab_t;
    typedef nt2::table<la_int,nt2::matlab_index_>   fitab_t;
    typedef nt2::table<type_t,index_t>                tab_t;
    typedef nt2::table<btype_t,index_t>              btab_t;
    typedef nt2::table<la_int,index_t>               itab_t;
      
      ////////////////////////////////////////////////////////////////////////////
      // General LU Solver
      //  A is            N x N
      //  X is            N x nrhs
      ////////////////////////////////////////////////////////////////////////////
      
      static inline bool LUSolveIP(XPR0&  A,
                                   XPR1&  X)
      {
        std::cout << nt2::type_id<itab_t>() << std::endl; 
        std::cout << "icitte1 " << std::endl;
        BOOST_ASSERT_MSG(nt2::issquare(A), "matrix A is not square");
        BOOST_ASSERT_MSG(nt2::areofsameheight(A, X), "A and X have different heights");
        std::cout << "icitte2 " << std::endl;
        long int info;
        const long int Ml   = nt2::height(A);
        const long int K    = nt2::width(X);
        const long int lda  = nt2::leading_size(A);
        const long int ldx  = nt2::leading_size(X);
        std::cout << "icitte3 " << Ml << " " << K << "  " << lda << "  " << ldx <<  "  " << std::endl;
        itab_t ipiv(nt2::of_size(Ml, 1)); 
        std::cout << "icitte4 " << std::endl;
        nt2::details::gesv (&Ml, &K, A.raw(), &lda, ipiv.raw(), X.raw(), &ldx, &info);
        std::cout << "icitte5 " << std::endl;
        std::cout << "info " << info << std::endl;
        
        for(int i=0; i < size(X, 1); i++)
          {
            std::cout << X(i) << ", "; 
          }
        std::cout << std::endl; 
        //        BOOST_ASSERT_MSG(info!= 0, "Lapack error : gesv in LUSolveIP");
        std::cout << "icitte6 " << std::endl;
        return (info == 0);
        
      }
      
      ////////////////////////////////////////////////////////////////////////////
      // General QR solver
      //  A is            M x N
      //  X is or will be N x nrhs
      //  B is            M x nrhs    
      ////////////////////////////////////////////////////////////////////////////
      static inline bool QRSolveIP(XPR0& A,
                                   XPR1& X,
                                   const XPR1& B)
      {

      // General QR solver
        std::cout << "icitte1 " << std::endl;
        BOOST_ASSERT_MSG(nt2::areofsameheight(A, B), "A and B have different heights");
        long int info;
        const long int Ml = size(A, 1);
        const long int Nl = size(A, 2);
        const int nrhs = size(B, 2);
        const long int nrhsl = nrhs;
        const long int lda = leading_size(A); 
        const char trans = 'n';
         std::cout << "icitte2 " << std::endl;
       
        // typically A is non-square, so we need to create tmp X because is
        //  X is N x nrhs, while B is M x nrhs.  We need to make copies of
        //  these so that the routine won't corrupt data around X and B
        
        if (Ml != Nl)
          {
            int Mm =  nt2::max(nt2::max(Ml,Nl),1l);
            table<type_t> Xtmp = B; //nt2::expand(B, nt2::of_size(Mm, nrhs));
            long int ldx = leading_size(Xtmp); 
            nt2::details::gels (&trans, &Ml, &Nl, &nrhsl,
                                A.raw(), &lda, Xtmp.raw(), &ldx, &info);
            X = Xtmp; //(Range(1, Nl), Range(1, nrhs)); 
            BOOST_ASSERT_MSG(info!= 0, "Lapack error : gels in QRSolveIP(1)");
          }
        else
          {
            std::cout << "icitte3 " << std::endl;
            X = B; 
            std::cout << "icitte4 " << std::endl;
            long int ldx = leading_size(X); 
            std::cout << "icitte5 " << std::endl;
            nt2::details::gels (&trans, &Ml, &Nl, &nrhsl,
                                A.raw(), &lda, X.raw(), &ldx, &info);
            std::cout << "icitte6 " << std::endl;
            BOOST_ASSERT_MSG(info!= 0, "Lapack error : gels in QRSolveIP(2)");
           }
        return info == 0; 
      }
      
      
      ////////////////////////////////////////////////////////////////////////////
      // General Cholevski solver
      //  A is            N x N
      //  X is or will be M x nrhs
      //  B is            M x nrhs
      //  need A symetric definite positive
      ////////////////////////////////////////////////////////////////////////////
      static inline bool CholSolveIP(arg0 A,
                                     arg1 X)
      {
        BOOST_ASSERT_MSG(nt2::issquare(A), "matrix A is not square");
        BOOST_ASSERT_MSG(nt2::areofsameheight(A, X), "A and X have different heights");
        long int info;
        const long int M = height(A);
        const long int K = width(X);
        const long int lda = leading_size(A);
        const long int ldx = leading_size(X);
        const char uplo = 'L';
        nt2::details::posv (&uplo, &M, &K,A.raw(), &lda, X.raw(), &ldx, &info);
        BOOST_ASSERT_MSG(info!= 0, "Lapack error : gels in CholSolveIP");
        return info == 0; 
      }
      
      
      ////////////////////////////////////////////////////////////////////////////
      // General SVD solver
      //  A is            M x N            may be rank-deficient
      //  X is or will be N x nrhs
      //  B is            M x nrhs    
      ////////////////////////////////////////////////////////////////////////////
      static inline long int SVDSolveIP(arg0 A,
                                        arg1 X,
                                        const arg2 B)
      {
        if(!areofsamedimensions(X, B)) { X.resize(size(B));   }
        long int info;
        const long int Ml = size(A, 1);
        const long int Nl = size(A, 2);
        const long int nrhs = size(X, 2);
        const long int lda = leading_size(A); 
        rtab_t S(of_size(nt2::min(Ml, Nl), 1)); 
        long int rank;
        rtype_t rcond =  Mone<rtype_t>();                                       
        
        // typically is A non-square, so we need to create tmp X because is 
        //  X is N x nrhs, while B is M x nrhs.  We need to make copies of
        //  these so that the routine won't corrupt data around X and B
        
        if (Ml != Nl)
          {
            long int Mm =  std::max(std::max(Ml,Nl),1l);
            tab_t Xtmp = B; //nt2::expand(B, nt2::of_size(Mm, nrhs));
            long int ldx = leading_size(Xtmp); 
            nt2::details::gelsd(&Ml, &Nl, &nrhs, A.raw(), &lda, Xtmp.raw(), &ldx,
                  S.raw(), &rcond, &rank, &info);                   
            X = Xtmp; //(Range(1, Nl), Range(1, nrhs));
            BOOST_ASSERT_MSG(info!= 0, "Lapack error : gelsd in SVDSolveIP(1)");
            return rank; 
          }
        else
          {
            X = B; 
            long int ldx = leading_size(X); 
            nt2::details::gelsd(&Ml, &Nl, &nrhs, A.raw(), &lda, X.raw(), &ldx,
                  S.raw(), &rcond, &rank, &info);                   
            BOOST_ASSERT_MSG(info!= 0, "Lapack error : gelsd in SVDSolveIP(2)");
            return rank; 
          }
      }
      
      ////////////////////////////////////////////////////////////////////////////
      //  Solve triangular systems
      //  A is            N x N
      //  X is            N x nrhs
      ////////////////////////////////////////////////////////////////////////////
      
      static inline bool TRISolveIP(arg0 A,
                                    arg1 X,
                                    const char &uplo,
                                    const char &t,
                                    const char &d )
      {
        BOOST_ASSERT_MSG(nt2::areofsameheight(A, X), "A and X have different heights");
        long int info;
        const long int Ml   = height(A);
        const long int K    = width(X);
        const long int lda  = leading_size(A); 
        const long int ldx  = leading_size(X); 
        nt2::details::trtrs((char*)&uplo, (char*)&t, (char*)&d, &Ml, &K, A.raw(), 
                            &lda, X.raw(), &ldx, &info);
        BOOST_ASSERT_MSG(info!= 0, "Lapack error : trtrs in TRISolveIP");
        return (info == 0); 
      }
    };
  }
}
#endif 
