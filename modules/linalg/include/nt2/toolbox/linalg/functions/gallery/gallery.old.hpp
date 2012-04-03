#ifndef __PROJECT__NT2__FILE__NT2_ALGEBRA_CORE_GALLERY_HPP__INCLUDED
#define __PROJECT__NT2__FILE__NT2_ALGEBRA_CORE_GALLERY_HPP__INCLUDED

////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2006 for LASMEA UMR 6602 du CNRS.                  
//  All rights reserved.                                             
//                                                                   
//  This file is part of the NT2 C++ Library.  This library is       
//  free software; you can redistribute it and/or modify it under    
//  the terms of the GNU Lesser General Public License as published  
//  by the Free Software Foundation; either version 2.1, or (at      
//  your option) any later version.                                  
//                                                                   
//  This library is distributed in the hope that it will be useful,  
//  but WITHOUT ANY WARRANTY; without even the implied warranty of   
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    
//  GNU Lesser General Public License for more details.              
//                                                                   
//  You should have received a copy of the GNU Lesser General        
//  Public License along with this library; see the file COPYING.    
//  If not, send mail to the developers of NT2                       
//                                                                   
//  As a special exception, you may use this file as part of a free  
//  software library without restriction.  Specifically, if other    
//  files instantiate templates or use macros or inline functions    
//  from this file, or you compile this file and link it with other  
//  files to produce an executable, this file does not by itself     
//  cause the resulting executable to be covered by the GNU Lesser   
//  General Public License.  This exception does not however         
//  invalidate any other reasons why the executable file might be    
//  covered by the GNU Lesser General Public License.                
//                                                                   
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//  what        : gallery.hpp                                        
//  who         : contributed by Joel FALCOU and Jean-Thierry LAPRESTE 
//  when        : Mon Sep 11 21:29:24 2006                           
//  where       : nt2/algebra/core/gallery.hpp                       
//  from        :                                                    
//  to          :                                                    
//  description :                                                    
//  modified    :                                                    
////////////////////////////////////////////////////////////////////////////////

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
    using  nc::matrix;

    //    Vandermonde matrix.  
    template < typename XPR > 
    inline matrix < typename ttt::float_promotion<typename XPR::type_t >::type_t >
    vander( const XPR & x1,  size_t w = size_t(-1))
    {
      typedef typename ttt::float_promotion<typename XPR::type_t >::type_t type_t;
      //      nc::matrix < double > x; 
      //      x = colvect(x1); 
      if (w ==  size_t(-1)) w =  x1.numel(); 
      return pow(colvect(x1)(All(), Repeat(Begin(), w)), fliplr(ci(x1.numel(),w)));
    }

    inline  matrix < ptrdiff_t > magic(const size_t &n){ 
    // %MAGIC  Magic square.
    // %   MAGIC(N) is an N-by-N matrix constructed from the integers
    // %   1 through N^2 with equal row, column, and diagonal sums.
    // %   Produces valid magic squares for all N > 0 except N = 2.
    
      if (n%2 == 1)   // % Odd order.
        {
          matrix < double > I, J; 
          meshgrid(iota(1, n), iota(1, n), J, I);
          matrix < double > A = mod(I+J-(n+3)/2,n);
          matrix < double > B = mod(I+2*J-2,n);
         return eve::convert < ptrdiff_t > (n*A + B + 1);
        }
      else if (n%4 == 0)   // % Doubly even order.
        {
          matrix < ptrdiff_t > I, J; 
          meshgrid(iota(1, n), iota(1, n), J, I);
          matrix < bool > K = fix(mod(I,4)/2) == fix(mod(J,4)/2);
          matrix < double > M = trans(reshape(iota(1, sSqr(n)),n,n));
          M(K) = sSqr(n) + 1 - M(K);
          return eve::convert < ptrdiff_t > (M); 
        }
      else  // % Singly even order.
        {
          size_t  p = n/2;
          matrix < double > M = magic(p);
          M = eve::catv( eve::cat(M)(M+2*sSqr(p))(),
                    eve::cat(M+3*sSqr(p))(M+sSqr(p))());
          if (n == 2) return  eve::convert < ptrdiff_t > (M); 
          matrix < ptrdiff_t > i = colvect(iota(1, p));
          size_t  k = (n-2)/4;
          matrix < ptrdiff_t > j = eve::cath(iota(1, k), iota((n-k+2), n));
          M(eve::catv(i, i+p),j) = M(eve::catv(i+p, i),j);
          i = eve::cons(k+1);
          j = eve::cat(1)(i)();
          M(eve::catv(i, i+p),j) = M(eve::catv(i+p, i),j);
          return  eve::convert < ptrdiff_t > (M);
        }
      return zeros(n); 
    }

    matrix < ptrdiff_t > inline rosser()
    {
      // %   The matrix is 8-by-8 with integer elements.
      // %   It has:
      // %       * A double eigenvalue.
      // %       * Three nearly equal eigenvalues.
      // %       * Dominant eigenvalues of opposite sign.
      // %       * A zero eigenvalue.
      // %       * A small, nonzero eigenvalue.
      ptrdiff_t r[] = {
        611,  196, -192,  407, -8,  -52,  -49,   29, 
        196,  899,  113, -192,  -71,  -43,   -8,  -44, 
        -192,  113,  899,  196, 61,   49,    8,   52, 
        407, -192,  196,  611,    8,   44,   59,  -23, 
        -8,  -71,   61,    8,  411, -599,  208,  208, 
        -52,  -43,   49,   44, -599,  411,  208,  208, 
        -49,  -8,    8,   59,  208,  208,   99, -911, 
        29,  -44,   52,  -23,   208,  208,  -911,   99};
      return matrix < ptrdiff_t > (&r[0], ofMatrixSize(8, 8)); 
      //        eve::cat(eve::cons(611.0,  196.0, -192.0,  407.0,   -8.0,  -52.0,  -49.0,   29.0))
      //         [eve::cons(196.0,  899.0,  113.0, -192.0,  -71.0,  -43.0,   -8.0,  -44.0)] 
      //         [eve::cons(-192.0,  113.0,  899.0,  196.0,   61.0,   49.0,    8.0,  52.0)]
      //         [eve::cons(407.0, -192.0,  196.0,  611.0,    8.0,   44.0,   59.0,  -23.0)]
      //         [eve::cons(-8.0,  -71.0,   61.0,    8.0,  411.0, -599.0,  208.0,  208.0)]
      //         [eve::cons(-52.0,  -43.0,   49.0,   44.0, -599.0,  411.0,  208.0,  208.0)]
      //         [eve::cons(-49.0,   -8.0,    8.0,   59.0,  208.0,  208.0,   99.0, -911.0)]
      //         [eve::cons(29.0,  -44.0,   52.0,  -23.0,  208.0,  208.0, -911.0,   99.0)]();
    }
    // /////////////////////////////////////////////////////////////////////////
    //     hilb   Hilbert matrix.
    // /////////////////////////////////////////////////////////////////////////
    inline  matrix < double > hilb(const size_t& n){
      return 1.0/(rif(n)+cif(n)-1);
    }
    
    // /////////////////////////////////////////////////////////////////////////
    //     cauchy    Cauchy matrix.
    // /////////////////////////////////////////////////////////////////////////
    template < class T > inline
    matrix < T > cauchy(const matrix < T > &x, const matrix < T > &y){
      //CAUCHY Cauchy matrix.
      //   c = cauchy(x,y), where X and Y are N-vectors, is the
      //   N-by-N matrix with C(i,j) = 1/(X(i)+Y(j)). By default, Y = X.
      //   If X is a scalar, GALLERY('CAUCHY',X) is the same as
      //   GALLERY('CAUCHY',1:X).
      //
      //   Explicit formulas are known for the elements of INV(C) and DET(C).
      //   DET(C) is nonzero if X and Y both have distinct elements.
      //   C is totally positive if 0 < X(1) < ... < X(N) and
      //                            0 < Y(1) < ... < Y(N).
      
      //   References:
      //   [1] D. E. Knuth, The Art of Computer Programming, Volume 1,
      //       Fundamental Algorithms, third edition, Addison-Wesley, Reading,
      //       Massachusetts, 1997.
      //   [2] E. E. Tyrtyshnikov, Cauchy-Toeplitz matrices and some applications,
      //       Linear Algebra and Appl., 149 (1991), pp. 1-18.
      //   [3] O. Taussky and M. Marcus, Eigenvalues of finite matrices, in
      //       Survey of Numerical Analysis, J. Todd, ed., McGraw-Hill, New York,
      //       1962, pp. 279-313. (The totally positive property is on p. 295.)
      //
      return T(1)/((ones(y.numel(), 1)*trans(y(All())))+(x(All())*ones(1, x.numel())));
    }
    
    //////////////////////////////////////////////////////////////////////////////
    template < class T >
    matrix < T > cauchy(const matrix < T > &x){
      return T(1)/((ones(size(x, 0), 1)* trans(x(All())))+(x(All())* ones(1, size(x, 0))));
      
    }
    
    //////////////////////////////////////////////////////////////////////////////
    //     chebspec  Chebyshev spectral differentiation matrix.
    template < class T > inline
    matrix < T > chebspec(size_t n, size_t k = 0){
      // CHEBSPEC Chebyshev spectral differentiation matrix.
      //    C = GALLERY('CHEBSPEC',N,K) is a Chebyshev spectral
      //       differentiation matrix of order N.  K = 0 (the default) or 1.
      //    For K = 0 ("no boundary conditions"), C is nilpotent, with
      //       C^N = 0 and it has the null vector ONES(N,1).
      //       C is similar to a Jordan block of size N with eigenvalue zero.
      //    For K = 1, C is nonsingular and well conditioned, and its
      //       eigenvalues have negative real parts.
      //    For both K, the computed eigenvector matrix X from EIG is
      //       ill-conditioned (MESH(REAL(X)) is interesting).
      
      //    References:
      //    [1] C. Canuto, M. Y. Hussaini, A. Quarteroni and T. A. Zang,
      //        Spectral Methods in Fluid Dynamics, Springer-Verlag, Berlin,
      //        1988, p. 69.
      //    [2] L. N. Trefethen and M. R. Trummer, An instability phenomenon in
      //        spectral methods, SIAM J. Numer. Anal., 24 (1987), pp. 1008-1023.
      //    [3] D. Funaro, Computing the inverse of the Chebyshev collocation
      //        derivative, SIAM J. Sci. Stat. Comput., 9 (1988), pp. 1050-1057.
      
      //    Nicholas J. Higham, Dec 1999.
      //    Copyright 1984-2002 The MathWorks, Inc. 
      //    $Revision: 1.10 $  $Date: 2002/04/15 03:41:53 $
      
      // k = 1 case obtained from k = 0 case with one bigger n.
      if(k == 1) n++;
      n--;
      matrix < T > c, d, one, x;
      c = zeros(n+1, n+1);
      one = ones(n+1,1);
      T k1 = pi<T>()/T(n); 
      x = cos(trans(iota<T>(0, 1, n))*k1);
      d = ones(n+1,1);
      d(0) = T(2); d(n) = T(2);
      // eye(size(C)) on next line avoids div by zero.
      c = div( ((d* trans(T(1)/d))), (x*trans(one))-(one*trans(x)) + eye(size(c)));
      //  Now fix diagonal and signs.
      c(0) = (2*n*n+1)/T(6);
      matrix < size_t > rr, cc;
      rr =  ri(c.height(), c.width())%2;
      cc =  ci(c.height(), c.width())%2; 
      //      c = where(x_or(rr, cc), -c, c); 
      for (size_t i=1; i <n;  i++){
        c(i,i) = -x(i)/(2.0*(1.0-std::pow(x(i), 2)));
      }
      c(n,n) = -c(0);
      if (k == 1){
        matrix < double > c1 =  c(Range(1, n), Range(1, n));
        return c1; 
      }
      return c; 
   }
    
    
    //////////////////////////////////////////////////////////////////////////////
    //     chebvand  Vandermonde-like matrix for the Chebyshev polynomials.
    template < class T > inline
    matrix < T > chebvand(size_t n){
      //        function C = chebvand(m,p)
      // CHEBVAND Vandermonde-like matrix for the Chebyshev polynomials.
      //    C = CHEBVAND(P), where P is a vector, produces the
      //    (primal) Chebyshev Vandermonde matrix based on the points P:
      //       C(i,j) = T_{i-1}(P(j)), where T_{i-1} is the Chebyshev
      //       polynomial of degree i-1.
      //    CHEBVAND(M,P) is a rectangular version of
      //    CHEBVAND(P) with M rows.
      //    Special case: If P is a scalar, then P equally spaced points on
      //       [0,1] are used.
      
      //    Reference:
      //    [1] N. J. Higham, Stability analysis of algorithms for solving confluent
      //        Vandermonde-like systems, SIAM J. Matrix Anal. Appl., 11 (1990),
      //        pp. 23-41.
      // 
      //    Nicholas J. Higham, Dec 1999.
      //    Copyright 1984-2002 The MathWorks, Inc. 
      //    $Revision: 1.11 $  $Date: 2002/04/15 03:41:56 $
      
      matrix < T > p = linspace(T(0),T(1),n);
      matrix < T > c = ones(n,n);
      if (n == 1) return c;
      c(1, All()) = p;
      //      Use Chebyshev polynomial recurrence.
      for (size_t i = 2; i < n; i++){
        c(i, All()) = 2*mul(p, c(i-1, All()))-c(i-2, All());
      }
      return c; 
    }
    //////////////////////////////////////////////////////////////////////////////
    template < class T > inline
    matrix < T > chebvand(const matrix < T >  p,  size_t m = 0){
      size_t n = p.numel(); 
      m = m ? m : n;
      matrix < T > c = ones(m,n);
      if (m == 1) return c; 
      c(1, All()) = p; 
      //      Use Chebyshev polynomial recurrence.
      for (size_t i = 2; i < m; i++){
        c(i, All()) = 2*mul(p,c(i-1, All()))-c(i-2, All());
      }
      return c; 
    }
    
    
    //////////////////////////////////////////////////////////////////////////////
    //toeplitz matrix
    template < class XPR1,  class XPR2>  inline
    matrix < typename XPR1::type_t  > toeplitz( const XPR1 &c,  XPR2 r)
    {
      //   TOEPLITZ Toeplitz matrix.
      //   TOEPLITZ(C,R) is a non-symmetric Toeplitz matrix having C as its
      //   first column and R as its first row.   
      //
      //   TOEPLITZ(R) is a symmetric (or Hermitian) Toeplitz matrix.
      //
      //   See also HANKEL.
      
      //   Revised 10-8-92, LS - code from A.K. Booer.
      //   Copyright 1984-2002 The MathWorks, Inc. 
      //   $Revision: 5.11 $  $Date: 2002/04/15 03:45:12 $
      
      // if nargin < 2,
      // c(1) = conj(c(1)); r = c; c = conj(c); // set up for Hermitian Toeplitz
      // else
      typedef typename XPR1::type_t T; 
      if (r.first_index() != c.first_index())
      {
        //     cerr << "toeplitz:DiagonalConflict: Column wins diagonal conflict.\n" <<  endl; 
      }
      r = reshape(r, ofSize(r.numel(), 1));// force column structure
      size_t p = r.length();
      size_t m = c.length();
      //    matrix < T >  x = catv(r(Range(End(), -1, Begin()+1)), c(All()));//build vector of user data
      matrix < T >  x = catv( 
           r(Range(End(), -1, Begin()+1)),
           c(All())); 
      matrix < ptrdiff_t > t0  = trans(ri(m, p) +  fliplr(ci(m, p)));
      t0.reshape(ofSize(1, m*p));
      matrix < T > t(ofSize(m, p));
      t(All()) =  x(t0); 
      matrix <T> res = trans(t);
      return res; // actual data
    }
    
    template < class XPR1>  inline
    matrix < typename XPR1::type_t > toeplitz(const XPR1 &r)
    {
      matrix < typename XPR1::type_t > c;
      c = conj(r);
      return  toeplitz<XPR1>(c, r); 
    }
    
    
    //////////////////////////////////////////////////////////////////////////////
    // chow      Chow matrix -- a singular Toeplitz lower Hessenberg matrix.
    template < class T >  inline
    matrix < T > chow(size_t n,  T alpha = T(1), T delta =  T(0)){
      // %CHOW Chow matrix (singular Toeplitz lower Hessenberg matrix).
      // %   A = GALLERY('CHOW',N,ALPHA,DELTA) returns A such that
      // %      A = H(ALPHA) + DELTA*EYE, where H(i,j) = ALPHA^(i-j+1).
      // %   H(ALPHA) has p = FLOOR(N/2) zero eigenvalues, the rest being
      // %   4*ALPHA*COS( k*PI/(N+2) )^2, k=1:N-p.
      // %   Defaults: ALPHA = 1, DELTA = 0.
      
      // %   References:
      // %   [1] T. S. Chow, A class of Hessenberg matrices with known eigenvalues
      // %       and inverses, SIAM Review, 11 (1969), pp. 391-395.
      // %   [2] G. Fairweather, On the eigenvalues and eigenvectors of a class of
      // %       Hessenberg matrices, SIAM Review, 13 (1971), pp. 220-221.
      // %   [3] I. Singh, G. Poole and T. Boullion, A class of Hessenberg matrices
      // %       with known pseudoinverse and Drazin inverse, Math. Comp., 29 (1975),
      // %       pp. 615-619.
      // %
      // %   Nicholas J. Higham, Dec 1999.
      // %   Copyright 1984-2002 The MathWorks, Inc. 
      // %   $Revision: 1.10 $  $Date: 2002/04/15 03:41:59 $
      
      // if nargin < 3, delta = 0; end
      // if nargin < 2, alpha = 1; end
      matrix < T >  z1, z2;
      z2 =  zeros(1, n); z2(0) = alpha;  z2(1) = 1; 
      z1 = pow(alpha, iota(1, 1, int(n))); 
      return toeplitz( z1, z2 ) + delta*eye(n, n);
    }
    
    //     circul    Circulant matrix.
    template < class T >  inline
    matrix < T > circul(matrix  < T > &v){
      //CIRCUL Circulant matrix.
      //   C = GALLERY('CIRCUL',V) is the circulant matrix whose first row is V.
      //   A circulant matrix has the property that each row is obtained
      //   from the previous one by cyclically permuting the entries one step
      //   forward. It is a special Toeplitz matrix in which the diagonals
      //   "wrap round". If V is a scalar, then C = GALLERY('CIRCUL',1:V).
      //
      //   The eigensystem of C (N-by-N) is known explicitly. If t is an Nth
      //   root of unity, then the inner product of V with W = [1 t t^2 ... t^(N-1)]
      //   is an eigenvalue of C, and W(N:-1:1) is an eigenvector of C.
      
      //   Reference:
      //   [1] P.J. Davis, Circulant Matrices, John Wiley, 1977.
      //
      //   Copyright 1984-2002 The MathWorks, Inc. 
      //   $Revision: 1.11 $  $Date: 2002/04/15 03:42:02 $
      
      int n = length(v);
      v = reshape(v, v.numel(), 1);   // Make sure v is a row vector.
      matrix < size_t > ii = iota(n, -1, 1);
      ii(0) = 0;
      matrix < T > v1 = v(ii); 
      return toeplitz( v1, v );
    }
    
    //     clement   Clement matrix -- tridiagonal with zero diagonal entries.
    template < class T >  inline
    matrix < T > clement(size_t n, size_t k = 0){
      //CLEMENT Clement matrix.
      //   GALLERY('CLEMENT',N,K) is a tridiagonal matrix with zero diagonal 
      //   entries and known eigenvalues. It is singular if N is odd. About 64
      //   percent of the entries of the inverse are zero. The eigenvalues
      //   are plus and minus the numbers N-1, N-3, N-5, ..., (1 or 0).
      //   For K = 0 (the default) the matrix is unsymmetric, while for
      //   K = 1 it is symmetric. GALLERY('CLEMENT',N,1) is diagonally similar 
      //   to GALLERY('CLEMENT',N).
      //
      //   Note:
      //   Similar properties hold for GALLERY('TRIDIAG',X,Y,Z) where 
      //   Y = ZEROS(N,1). The eigenvalues still come in plus/minus pairs but 
      //   they are not known explicitly.
      //
      
      //   References:
      //   [1] P.A. Clement, A class of triple-diagonal matrices for test
      //   purposes, SIAM Review, 1 (1959), pp. 50-52.
      //   [2] O. Taussky and J. Todd, Another look at a matrix of Mark Kac,
      //   Linear Algebra and Appl., 150 (1991), pp. 341-360.
      //
      //   Copyright 1984-2002 The MathWorks, Inc. 
      //   $Revision: 1.9 $  $Date: 2002/04/15 03:41:17 $
      
      n--;
      matrix < T >  t;
      matrix < size_t > x= iota(int(n), -1, 1);
      matrix < size_t > z =iota(1, 1, int(n));
      
      if (k == 0){
        t = diag(x, -1) + diag(z, 1);
      } else {
        matrix < T > y= sqrt(mul(x,z));
        t = diag(y, -1) + diag(y, 1);
      }; 
      return t;
    }
    
    //////////////////////////////////////////////////////////////////////////////
    //     condex    Counter-examples to matrix condition number estimators.
    // todo
    
    //////////////////////////////////////////////////////////////////////////////
    template < class T >  inline
    matrix < T > cycol(size_t m, size_t n, size_t k = 0){
      // function A = cycol(n, k)
      // CYCOL  Matrix whose columns repeat cyclically.
      //    A = GALLERY('CYCOL',[M N], K) is an M-by-N matrix of the form 
      //    A = B(1:M,1:N) where B = [C C C...] and C = RANDN(M, K). Thus A's 
      //    columns repeat cyclically, and A has rank at most K. K need not 
      //    divide N. K defaults to ROUND(N/4).
      //    GALLERY('CYCOL',N,K), where N is a scalar, is the same as 
      //    GALLERY('CYCOL',[N N], K).
      
      //    Note:
      //    This type of matrix can lead to underflow problems for Gaussian
      //    elimination. See the reference below.
      // 
      //    Reference:
      //    [1] NA Digest Volume 89, Issue 3 (January 22, 1989).
      // 
      //    Copyright 1984-2002 The MathWorks, Inc. 
      //    $Revision: 1.10 $  $Date: 2002/04/15 03:41:26 $
      
      // Parameter n specifies dimension: m-by-n.
      
      if(k == 0) k = size_t(std::max(sRound (n/4.0),1.0));
      matrix < T > a = randn(m, k); 
      matrix < T > c(ofSize(m, 0)),  c1; 
      size_t deb = 0;
      for (size_t i=1;  i <= sCeil(1.0*n/k);  i++){
        c1 = cath(c, a(All(),AllTo(k-1)));
        c =  c1; 
        deb+= k; 
      }
     a = c(All(), AllTo(n-1));
      return a; 
    }
    
    
    
    //     tridiag   Tridiagonal matrix (sparse).
    template < class XPR1,  class XPR2,  class XPR3>  inline
    matrix < double > tridiag(const XPR1 & x, const XPR2 & y, const XPR3 & z){
      matrix < double > res = diag(x(All()), -1) + diag(y(All())) + diag(z(All()), 1); 
      return res; 
    }
    
    
    
    
    
    //     dorr      Dorr matrix -- diagonally dominant, ill-conditioned, tridiagonal.
    //               (One or three output arguments, sparse)
    template < class T >  inline
    matrix < T > dorr(size_t n, T theta = 0.01){
      // DORR Dorr matrix (sparse).
      //    [C,D,E] = GALLERY('DORR',N,THETA) returns the vectors defining a row
      //    diagonally dominant, tridiagonal N-by-N matrix that is ill-conditioned
      //    for small values of THETA >= 0. THETA defaults to 0.01.
      // 
      //    A = GALLERY('DORR',N,THETA) returns the (sparse) Dorr matrix itself.  
      //    This is the same as
      //        [C,D,E] = GALLERY('DORR',N,THETA);
      //        A = GALLERY('TRIDIAG',C,D,E);
      // 
      //    The columns of INV(C) vary greatly in norm.
      //    The amount of diagonal dominance, ignoring rounding errors, is:
      //         GALLERY('COMPAR',C)*ONES(N,1) = THETA*(N+1)^2 * [1 0 0 ... 0 1]'.
      
      //    Reference:
      //    [1] F. W. Dorr, An example of ill-conditioning in the numerical
      //        solution of singular perturbation problems, Math. Comp.,
      //        25 (1971), pp. 271-283.
      // 
      //    Nicholas J. Higham, Dec 1999.
      //    Copyright 1984-2002 The MathWorks, Inc. 
      //    $Revision: 1.11 $  $Date: 2002/04/15 03:41:29 $
      
      matrix < double, settings < FORTRAN_policy > > c =  zeros(n, 1),  d, e;
      d = c;
      e = c; 
      // c = zeros(n,1); e = c; d = c;
      // % All length n for convenience.  Make c, e of length n-1 later.
      T h = 1.0/(n+1);
      size_t m = size_t(sFloor ( (n+1)/2.0 ));
      T term = theta/sSqr(h);

      matrix < size_t > i = iota(1, m);
      c(i) = -term;
      e(i) = c(i) - (0.5-i*h)/h;
      d(i) = -(c(i) + e(i));
      
      i = iota(m+1, n);
      e(i) = -term;
     c(i) = e(i) + (0.5-i*h)/h;
      d(i) = -(c(i) + e(i));
      
      matrix < T > c1 = c(AllFrom(2));
      matrix < T > e1 = e(AllTo(n-1));
      return  tridiag (c1, d, e1);
      
    }
    
    // //     dramadah  Matrix of ones and zeroes whose inverse has large integer entries.
    template < class T >  inline
    matrix < T > dramadah(int n, size_t k = 1)
    {
      // function A = dramadah(n, k)
      // %DRAMADAH Matrix of zeros and ones with large determinant or inverse.
      // %   A = GALLERY('DRAMADAH',N,K) is an N-by-N (0,1) matrix for which
      // %   MU(A) = NORM(INV(A),'FRO') or DET(A) is relatively large.
      // %
      // %   K = 1: (default)
      // %      A is Toeplitz, with ABS(DET(A)) = 1, and MU(A) > c(1.75)^N,
      // %      where c is a constant. INV(A) has integer entries.
      // %   K = 2:
      // %      A is upper triangular and Toeplitz. INV(A) has integer entries.
      // %   K = 3:
      // %      A has maximal determinant among (0,1) lower Hessenberg matrices.
      // %      DET(A) = the n'th Fibonacci number. A is Toeplitz.
      // %      The eigenvalues have an interesting distribution in the complex
      // %      plane.
      // %
      // %   An anti-Hadamard matrix A is a matrix with elements 0 or 1 for
      // %   which MU(A) = NORM(INV(A),'FRO') is maximal.  For K = 1,2 this function
      // %   returns matrices with MU(A) relatively large, though not necessarily
      // %   maximal.
      
      // %   References:
      // %   [1] R. L. Graham and N. J. A. Sloane, Anti-Hadamard matrices,
      // %       Linear Algebra and Appl., 62 (1984), pp. 113-137.
      // %   [2] L. Ching, The maximum determinant of an nxn lower Hessenberg
      // %       (0,1) matrix, Linear Algebra and Appl., 183 (1993), pp. 147-153.
      // %
      // %   Nicholas J. Higham, Dec 1999.
      // %   Copyright 1984-2002 The MathWorks, Inc. 
      // %   $Revision: 1.11 $  $Date: 2002/04/15 03:42:08 $
      
      if(k == 1)
      {
        // Toeplitz
        matrix < T, settings<FORTRAN_policy> > c = ones(n,1);
        for(int i=2; i <= n; i+= 4)
        {
          size_t m = std::min(1,n-i);
          c(Range (i, i+m)) = T(0);
        }
        matrix < T,settings<FORTRAN_policy> > r = zeros(n,1);
        r(1) = r(2) = r(4) = 1; 
        if (n < 4) r = r(Range (1, n));
        return toeplitz (c,r);
      } 
      else if (k == 2)
      { // Upper triangular and Toeplitz
        matrix < T, settings<FORTRAN_policy> > c = zeros(n,1);
        c(1) = 1;
        matrix < T, settings<FORTRAN_policy> > r = ones(n,1);
        r(Range (3, 2, n)) = T(0);
        return toeplitz (c,r);
      } 
      else 
      { //   Lower Hessenberg.
        matrix < T, settings<FORTRAN_policy> > c = ones(n,1);
        c(Range (2, 2, n)) = T(0);
        matrix < T, settings<FORTRAN_policy> > d = zeros(n,1);
        d(1) = d(2) = 1; 
        return toeplitz  (c, d);
      }
      
    }
    
    
    
    //   //     fiedler   Fiedler matrix -- symmetric.
    // function A = fiedler(c)
    template < class XPR >
    matrix <typename XPR::type_t > fiedler(const XPR & c){
      // %FIEDLER Fiedler matrix.
      // %   A = GALLERY('FIEDLER',C), where C is an N-vector, is the N-by-N 
      // %   symmetric matrix with elements ABS(C(i)-C(j)). If C is a scalar, 
      // %   then A = GALLERY('FIEDLER',1:C).
      // %
      // %   A has a dominant positive eigenvalue and all the other eigenvalues
      // %   are negative. (Szego 1936)
      // %
      // %   Note: Explicit formulas for INV(A) and DET(A) are given in (Todd 1977)
      // %   and attributed to Fiedler. These indicate that INV(A) is 
      // %   tridiagonal except for nonzero (1,n) and (n,1) elements.
      
      // %   References:
      // %     [1] G. Szego, Solution to problem 3705, Amer. Math. Monthly,
      // %       43 (1936), pp. 246-259.
      // %     [2] J. Todd, Basic Numerical Mathematics, Vol. 2: Numerical Algebra,
      // %       Birkhauser, Basel, and Academic Press, New York, 1977, p. 159.
      // %
      // %   Copyright 1984-2002 The MathWorks, Inc. 
      // %   $Revision: 1.11 $  $Date: 2002/04/15 03:42:11 $
      
      size_t n = numel(c);
      
      // %  Handle scalar c.
      //  if (n == 1)
      //  n = c;
      //    c = 1:n;
      // end
      typedef typename XPR::type_t type_t; 
      matrix < type_t > d = reshape(c, 1, c.numel());//                    % Ensure c is a row vector.
      //    matrix < type_t > a  = abs(d(zeros(n,1), All())-trans(d(zeros(n,1), All())));
      matrix < type_t > a  = abs(d(Repeat(Begin(), n),  All())-trans(d(Repeat(Begin(), n), All())));
      return a; 
    }

 
    template < class T> 
    matrix < T > wilkinson(const size_t & n)
    {
      // %WILKINSON Wilkinson's eigenvalue test matrix.
      // %   WILKINSON(n) is J. H. Wilkinson's eigenvalue test matrix, Wn+.
      // %   It is a symmetric, tridiagonal matrix with pairs of nearly,
      // %   but not exactly, equal eigenvalues.
      // %   The most frequently used case is WILKINSON(21).
      // %   For example, WILKINSON(7) is
      // %
      // %          3  1  0  0  0  0  0
      // %          1  2  1  0  0  0  0
      // %          0  1  1  1  0  0  0
      // %          0  0  1  0  1  0  0
      // %          0  0  0  1  1  1  0
      // %          0  0  0  0  1  2  1
      // %          0  0  0  0  0  1  3
      // %
      T m = (n-1)/2;
      matrix < T > e = ones<T>(n-1,1);
      return diag(abs(iota(-m, m))) + diag(e,1) + diag(e,-1);
    }
    
    //     forsythe  Forsythe matrix -- a perturbed Jordan block.
    //     frank     Frank matrix -- ill-conditioned eigenvalues.
    //     gearmat   Gear matrix.
    //     grcar     Grcar matrix -- a Toeplitz matrix with sensitive eigenvalues.
    //     hanowa    Matrix whose eigenvalues lie on a vertical line in the complex
    //               plane.
    //     house     Householder matrix. (Three output arguments)
    //     invhess   Inverse of an upper Hessenberg matrix.
    //     invol     Involutory matrix.
    //     ipjfact   Hankel matrix with factorial elements. (Two output arguments)
    //     jordbloc  Jordan block matrix.
    //     kahan     Kahan matrix -- upper trapezoidal.
    //     kms       Kac-Murdock-Szego Toeplitz matrix.
    //     krylov    Krylov matrix.
    //     lauchli   Lauchli matrix -- rectangular.
    //     lehmer    Lehmer matrix -- symmetric positive definite.
    //     leslie    Leslie matrix.
    //     lesp      Tridiagonal matrix with real, sensitive eigenvalues.
    //     lotkin    Lotkin matrix.
    //     minij     Symmetric positive definite matrix MIN(i,j).
    //     moler     Moler matrix -- symmetric positive definite.
    //     neumann   Singular matrix from the discrete Neumann problem (sparse).
    //     orthog    Orthogonal and nearly orthogonal matrices.
    //     parter    Parter matrix -- a Toeplitz matrix with singular values near PI.
    //     pei       Pei matrix.
    //     poisson   Block tridiagonal matrix from Poisson's equation (sparse).
    //     prolate   Prolate matrix -- symmetric, ill-conditioned Toeplitz matrix.
    //     randcolu  Random matrix with normalized cols and specified singular values.
    //     randcorr  Random correlation matrix with specified eigenvalues.
    //     randhess  Random, orthogonal upper Hessenberg matrix.
    //     randjorth Random J-orthogonal matrix.
    //     rando     Random matrix with elements -1, 0 or 1.
    //     randsvd   Random matrix with pre-assigned singular values and specified
    //               bandwidth.
    //     redheff   Matrix of 0s and 1s of Redheffer.
    //     riemann   Matrix associated with the Riemann hypothesis.
    //     ris       Ris matrix -- a symmetric Hankel matrix.
    //     smoke     Smoke matrix -- complex, with a "smoke ring" pseudospectrum.
    //     toeppd    Symmetric positive definite Toeplitz matrix.
    //     toeppen   Pentadiagonal Toeplitz matrix (sparse).
    
    //     triw      Upper triangular matrix discussed by Wilkinson and others.
    //     wathen    Wathen matrix -- a finite element matrix (sparse, random entries).
    //     wilk      Various specific matrices devised/discussed by Wilkinson.
    //               (Two output arguments)
    
    //     GALLERY(3) is a badly conditioned 3-by-3 matrix.
    //     GALLERY(5) is an interesting eigenvalue problem.  Try to find
    //     its EXACT eigenvalues and eigenvectors.
  }

  // ///////////////////////////////////////////////////////////////////////////
  //  End of alg namespace
  // ///////////////////////////////////////////////////////////////////////////

}

// /////////////////////////////////////////////////////////////////////////////
//  End of nt2 namespace
// /////////////////////////////////////////////////////////////////////////////


#endif

// /////////////////////////////////////////////////////////////////////////////
// End of gallery.hpp
// /////////////////////////////////////////////////////////////////////////////
