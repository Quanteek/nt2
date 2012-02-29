#error "This file is for documentation purpose only."
/*!
 * \ingroup toolbox
 * \defgroup algebra algebra functions
 * \brief Defines algebra functions
 **/
/*!
 * \ingroup algebra
 * \defgroup algebra
 * \brief Defines algebra
 **/
/*
  The algebra interface is based on 6 levels of interface from external
  libraries according to the language and the way those are written down, .

  Today blas/lapack 
  In the future other interesting external libraries as scalapack pblas linpack etc.

  1) declaration of C like calls translating the direct FORTRAN Call
  
  2) overloading of the C call in C++ over the elements type supportedl namely
     float,  double,  complex<float>  and complex<double>
     
  3) suppressing the shape parameters and remplacing array pointers by containers
     aware of their access properties
     grouping the interpretation parameters in a static (or not) type
     that will replace all of them (side uplo diag transa transb etc.) 

     Perhaps this level is to be divided in two parts

     Note that this can imply that our calls have more than 6 parameters

  4) An nt2 interface dealing with the input/output character of many fortran
     parameters (including the aliasing problem)

  5) Matlab like direct calls with no verification of adequation of the
     chosen routine to the matrix type/shape under the user responsability

  6) drivers interfaces allowing to use/select autmatically best adapted routines at
     run-time (as  mldivide in Matlab(tm) for instance)
     or syntax catching to replace patterns with unique call
     (for instance a*trans(A)*B+b*C by the proper gemm call)
     
  The two first levels are easy and quite automatic from what was written in lpp

  The third level depends on what can be obtained from the containers
          the shape of the matrix
          the fact that a container has to be considered transposed or not
          on its access
          the character symetric/hermitian triangular diag (meaning unit on the diag)
          to overide the storage

  The fourth level has to be driven somehow by the constness of the parameters

  The five level is quite straightforward and sugar to ease conversion matlab->nt2

  The sixth level is linked to the knowledge of matlab selection algorithms in order
  to properly mimick its choices.

 ===============================================================================================
 Example from blas
 -----------------------------------------------------------------------------------------------
 level 1 C interface: 30 calls (with c, z, s, d and ge sy he)
 defining 30 different calls as:
 
            extern "C"
              {
                void NT2_F77NAME(sgemm)(const char *transa, const char *transb, const long int *m, 
                      const long int *n, const long int *k, const float *alpha, 
                      const float *a, const long int *lda, const float *b, 
                      const long int *ldb, const float *beta, float *c, 
                      const long int *ldc);
              }
 -----------------------------------------------------------------------------------------------
 -----------------------------------------------------------------------------------------------
 level 2 C++ interface: 3 routines (with ge sy he) defining 10 functions as:
 
 #define NT2_MM(T, PREFIX)                                              \
 inline void gemm(const char *ta, const char *tb, const long int *m,    \
      const long int *n, const long int *k,                                  \
      const T *al, const T *a,                                               \
      const long int *lda, const T *b,                                       \
      const long int *ldb, const T *be, T *c,                                \
      const long int *ldc)                                                   \
 {                                                                      \
      BOOST_PP_CAT(PREFIX,BOOST_PP_CAT(gemm,_))(ta,tb,m,n,k,al,a,lda,b,ldb,be,c,ldc); \
 }                                                                      \
 
 NT2_MM(double, d) 
 NT2_MM(float,  s) 
 NT2_MM(std::complex<double>, z) 
 NT2_MM(std::complex<float>, c) 
 
 #undef NT2_MM
 -----------------------------------------------------------------------------------------------
 -----------------------------------------------------------------------------------------------
 level 3 C++ interface 1 common call for all routines  
 * namespace nt2
 * {
 *   template <class A0,class A1,class A2,class A3,class A4>
 *   void b_mm(const gem_status& gs,
 *             A0& c, A1 const& a, const A2 & b,
 *             A3 const& alpha = 1, A4 const& beta = 0);
 * }
 -----------------------------------------------------------------------------------------------
 -----------------------------------------------------------------------------------------------
 level 4 C++ interface 2 call for internal alias decision  
 * namespace nt2
 * {
 *   template <class A0,class A1,class A2,class A3,class A4>
 *   A0 b_mm(const gem_status& gs,
 *           const A0& c, A1 const& a, const A2 & b,
 *           A3 const& alpha = 1, A4 const& beta = 0);
 *
 *   template <class A0,class A1,class A2,class A3,class A4>
 *   A0 b_mm(const gem_status& gs,
 *           A0& c, A1 const& a, const A2 & b,
 *           A3 const& alpha = 1, A4 const& beta = 0);
 * }
 -----------------------------------------------------------------------------------------------
 -----------------------------------------------------------------------------------------------
 level 5 C++ interface is not relevant today for matrix multiplication (sparse matrices can
 further change this fact)
 -----------------------------------------------------------------------------------------------
 -----------------------------------------------------------------------------------------------
 level 6 C++ interface is to transform an expression of the form  D = a*trans(A)*B+b*C  to
 the proper b_mm call
 -----------------------------------------------------------------------------------------------
 -----------------------------------------------------------------------------------------------
 
 ===============================================================================================
 LAPACK shapes

 * These are the type/form of matrices supported by lapack in alphabetical order
 *
 * BD   bidiagonal
 * DI   diagonal
 * GB   general band
 * GE   general (i.e., unsymmetric, in some cases rectangular)
 * GG   general matrices, generalized problem (i.e., a pair of general matrices)
 * GT   general tridiagonal
 * HB   (complex) Hermitian band
 * HE   (complex) Hermitian
 * HG   upper Hessenberg matrix, generalized problem (i.e a Hessenberg and a triangular matrix)
 * HP   (complex) Hermitian, packed storage
 * HS   upper Hessenberg
 * OP   (real) orthogonal, packed storage
 * OR   (real) orthogonal
 * PB   symmetric or Hermitian positive definite band
 * PO   symmetric or Hermitian positive definite
 * PP   symmetric or Hermitian positive definite, packed storage
 * PT   symmetric or Hermitian positive definite tridiagonal
 * SB   (real) symmetric band
 * SP   symmetric, packed storage
 * ST   (real) symmetric tridiagonal
 * SY   symmetric
 * TB   triangular band
 * TG   triangular matrices, generalized problem (i.e., a pair of triangular matrices)
 * TP   triangular, packed storage
 * TR   triangular (or in some cases quasi-triangular)
 * TZ   trapezoidal
 * UN   (complex) unitary
 * UP   (complex) unitary, packed storage

 Ad hoc types must be designed to hold the  bidiagonal and tridiagonal matrices (symetric hermitian or not)
 as FORTRAN routines on bi-diag or tri-diag chose to pass several vectors instead of one array 

 Also not all prefixes deal with types GG for instance

 In Lapack (exerpt from the lapack doc)
 ----------
   An unsymmetric tridiagonal matrix of order n is stored in three one-dimensional arrays,
   one of length n containing the diagonal elements,
   and two of length n-1 containing the subdiagonal and superdiagonal elements in elements 1:n-1.

  A symmetric tridiagonal or bidiagonal matrix is stored in two one-dimensional arrays,
  one of length n containing the diagonal elements,
  and one of length n-1 containing the off-diagonal elements.
  (EISPACK routines store the off-diagonal elements in elements 2:n of a vector of length n.)
 ----------

 This must be replaced by some supplementary shapes 
