#ifndef __PROJECT__FILE__NT2_ALGEBRA_SOLVERS_INTERNAL_SOLVERS_HPP__INCLUDED
#define __PROJECT__FILE__NT2_ALGEBRA_SOLVERS_INTERNAL_SOLVERS_HPP__INCLUDED

////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2007 for LASMEA UMR 6602 du CNRS.                  
//  All rights reserved.                                             
//                                                                   
//  License informations are available in the LICENSE file.          
//  Additionnal informations are available on http://nt2.sourceforge.net 
////////////////////////////////////////////////////////////////////////////////
#include <lpp/lapack_cpp_itf/gesv_itf.hh>
#include <lpp/lapack_cpp_itf/posv_itf.hh>
#include <algorithm>


// /////////////////////////////////////////////////////////////////////////////
//  Beginning of nt2 namespace
// /////////////////////////////////////////////////////////////////////////////
namespace nt2
{
  
  // ///////////////////////////////////////////////////////////////////////////
  //  Beginning of alg namespace
  // ///////////////////////////////////////////////////////////////////////////
  namespace alg
  {
    using namespace lpp;
    
    template < class XPR0,  class XPR1,  class XPR2 = XPR1> 
    class lasolvers
    {
    public:
      typedef typename base_of<XPR0>::type                          base0;
      typedef typename base_of<XPR1>::type                          base1;
      typedef typename base_of<XPR2>::type                          base2;

      typedef typename base0::type                                  type0;
      typedef typename base1::type                                  type1;
      typedef typename base2::type                                  type2;
      
      typedef typename result<type0,type1,type2>::floating       floating;
      typedef typename result<type0,type1,type2>::real_float realfloating;

      typedef typename ref_from<XPR0>::type                          arg0;
      typedef typename ref_from<XPR1>::type                          arg1;
      typedef typename ref_from<XPR2>::type                          arg2;
      
      typedef Matricial_Checks<>                                   mcheck; 
      ////////////////////////////////////////////////////////////////////////////
      // General LU Solver
      //  A is            N x N
      //  X is            N x nrhs
      ////////////////////////////////////////////////////////////////////////////
      
      static inline bool LUSolveIP(arg0  A,
                                   arg1  X)
      {
        mcheck::SquareTest(__FILE__, __LINE__, A);
        mcheck::HeightsTest(__FILE__, __LINE__, A, X);
        long int info;
        const long int Ml   = A.height();
        const long int K    = X.width();
        const long int lda  = A.height();
        const long int ldx  = X.height();
        matrix < long int > ipiv(ofSize(Ml, 1)); 
        gesv (&Ml, &K, A.begin(), &lda, ipiv.begin(), X.begin(), &ldx, &info);
        mcheck::LapackTest(__FILE__, __LINE__, "gesv", A, info); 
        return (info == 0);
        
      }
      
      ////////////////////////////////////////////////////////////////////////////
      // General QR solver
      //  A is            M x N
      //  X is or will be N x nrhs
      //  B is            M x nrhs    
      ////////////////////////////////////////////////////////////////////////////
      static inline bool QRSolveIP(arg0 A,
                                   arg1 X,
                                   const arg2 B)
      {
        mcheck::HeightsTest(__FILE__, __LINE__, A, B);
        //        if(!have_same_dimensions(X, B)) { X.resize(B.size());   }
        long int info;
        const long int Ml = A.size(1);
        const long int Nl = A.size(2);
        const int nrhs = X.size(2);
        const long int nrhsl = nrhs;
        const long int lda = A.size(1); 
        const char trans = 'n';
        
        // typically A is non-square, so we need to create tmp X because is
        //  X is N x nrhs, while B is M x nrhs.  We need to make copies of
        //  these so that the routine won't corrupt data around X and B
        
        if (Ml != Nl)
          {
            int Mm =  std::max(std::max(Ml,Nl),1l);
            matrix < typename XPR2::type_t > Xtmp = expand(B, ofSize(Mm, nrhs));
            long int ldx = Xtmp.size(1); 
            gels (&trans, &Ml, &Nl, &nrhsl, A.begin(), &lda, Xtmp.begin(), &ldx, &info);
            X = Xtmp(Range(0, Nl-1), Range(0, nrhs-1)); 
          }
        else
          {
            X = B; 
            long int ldx = X.size(1); 
            gels (&trans, &Ml, &Nl, &nrhsl, A.begin(), &lda, X.begin(), &ldx, &info);
          }
        mcheck::LapackTest(__FILE__, __LINE__, "gels", A, info); 
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
        mcheck::SquareTest(__FILE__, __LINE__, A);
        mcheck::HeightsTest(__FILE__, __LINE__, A, X);
        //        if(!is_sym(A)){ return false; }
        long int info;
        const long int M = A.height();
        const long int K = X.width();
        const long int lda = A.height();
        const long int ldx = X.height();
        const char uplo = 'L';
        posv (&uplo, &M, &K,A.begin(), &lda, X.begin(), &ldx, &info);
        mcheck::LapackTest(__FILE__, __LINE__, "posv", A, info); 
        //               if (info > 0)std::cout << "A is not symmetric-positive-definite." << std::endl; 
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
        mcheck::HeightsTest(__FILE__, __LINE__, A, B);
        if(!have_same_dimensions(X, B)) { X.resize(B.size());   }
        long int info;
        const long int Ml = A.size(1);
        const long int Nl = A.size(2);
        const long int nrhs = X.size(2);
        const long int lda = A.size(1); 
        matrix <realfloating > S(ofMatrixSize(sMin(Ml, Nl), 1)); 
        long int rank;
        realfloating rcond =  -1.0;                                       
        
        // typically is A non-square, so we need to create tmp X because is 
        //  X is N x nrhs, while B is M x nrhs.  We need to make copies of
        //  these so that the routine won't corrupt data around X and B
        
        if (Ml != Nl)
          {
            long int Mm =  std::max(std::max(Ml,Nl),1l);
            matrix < typename XPR2::type_t > Xtmp = expand(B, ofSize(Mm, nrhs));
            long int ldx = Xtmp.size(1); 
            gelsd(&Ml, &Nl, &nrhs, A.begin(), &lda, Xtmp.begin(), &ldx,
                  S.begin(), &rcond, &rank, &info);                   
            X = Xtmp(Range(0, Nl-1), Range(0, nrhs-1));
            mcheck::LapackTest(__FILE__, __LINE__, "gelsd", A, info); 
            return rank; 
          }
        else
          {
            X = B; 
            long int ldx = X.size(1); 
            gelsd(&Ml, &Nl, &nrhs, A.begin(), &lda, X.begin(), &ldx,
                  S.begin(), &rcond, &rank, &info);                   
            mcheck::LapackTest(__FILE__, __LINE__, "gelsd", A, info); 
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
        mcheck::HeightsTest(__FILE__, __LINE__, A, X);
        long int info;
        const long int Ml   = A.height();
        const long int K    = X.width();
        const long int lda  = A.height();
        const long int ldx  = X.height();
        trtrs ((char*)&uplo, (char*)&t, (char*)&d, &Ml, &K, A.begin(), &lda, X.begin(), &ldx, &info);
        //        dtrsv_((char*)&uplo, (char*)&t, (char*)&d, &Ml, A.begin(), &lda, X.begin(), &ldx);
        mcheck::LapackTest(__FILE__, __LINE__, "trtrs", A, info); 
        return (info == 0); 
      }
            
    };
    
  }
  
  // ///////////////////////////////////////////////////////////////////////////
  //  End of alg namespace
  // ///////////////////////////////////////////////////////////////////////////
 
}

// /////////////////////////////////////////////////////////////////////////////
//  End of nt2 namespace
// /////////////////////////////////////////////////////////////////////////////


#endif // __PROJECT__FILE__NT2_ALGEBRA_SOLVERS_INTERNAL_SOLVERS_HPP__INCLUDED

// /////////////////////////////////////////////////////////////////////////////
// End of internal_solvers.hpp
// /////////////////////////////////////////////////////////////////////////////
