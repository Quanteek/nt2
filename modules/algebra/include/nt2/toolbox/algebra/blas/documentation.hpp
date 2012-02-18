#error "This file is for documentation purpose only."
/*!
 * \ingroup toolbox
 * \defgroup blas blas functions
 * \brief Defines blas functions
 **/
/*!
 * \ingroup blas
 * \defgroup blas
 * \brief Defines blas
 *
 * The aim of this part is to encapsulate the blas functions
 * in a generic manner suited to C++
 *
 * rationale
 * blas functions consits of name and parameters
 *
 * Blas names:
 * blas function have names consisting in 3 parts:
 * -- one charaxter among s, d, c, z defining the type of elements datas
 *    with which the routine will act namely float,  double,  complex<float>  and complex<double>
 *
 * -- one consisting of two characters dealing with the shape of involved matrices
 *    ge,  gb, 
 *    sy,  sb, sp,
 *    he,  hb, hp
 *    tr,  tb, tp
 *    In fact g, b, p correspond to a storage scheme: general, band, packed
 *            s,  h, t to a property of the matrix  :symetric, hermitian, triangular
 * The storage define how the elemnt are accessed in memory
 * The property tell what are the elements not accessed: except the g case,
 * almost half of the elements are unnessary to completely define the matrix
 *
 * -- finally one part devoted to define the kind of performed operation
 *    for instance mm means matrix matrix multiplication and mv means matrix vector
 *    multiplication
 *
 * Blas parameters
 * they are all addresses and of three species
 * -- chars that add information on the way to compute and on characteristics of the matrices
 *    These can be associated to the following properties:
 * .. Which are the pertinent values:
 *    * SIDE meaning the operation is a Left or Right operation
 *    * UPLO meaning the pertinent matix part is Upper or Lower
 *    * DIAG meaning the diagonal element are Unit or Not necesseraly unit
 * .. In which way the matrix must be considered:
 *    * trans mmeaning the matrix is Transposed Not transpose or Conjugated
 *
 * -- floating pointers to a matrix, vector or scalar data consisting of contiguous elements
 *    stored in the schemes defined by g, b, or p
 *
 * -- long int pointers containing matrices dimensions or extent
 *
 * Our routines will have name consisting of the two characters defining the kind
 * of performed operation and will use much less parameters than the original ones
 *
 * * The data pointers will be nt2 containers which are aware of their shape and of their
 *   dimensions and extents,  suppressing the need of passing dimensions as supplementary long ints
 *   parameters as well as storage information in the function name.
 *
 * * The char properties can be subsumed in a struct containing only static constants allowing to
 *   know the properties at compile time (but disallowing to chose to change dynamically the property, 
 *   we think this is rarely needed and can be replaced by testing and branching if not avoidable)
 *
 * So for example,  the functions gemm,  gbmm, hemm, hbmm, symm, sbmm
 * which all compute in various cases
 *                                C <- alpha*A°*B°+beta*C
 * with FORTRAN definition call analog to the C
 *    void gemm(const char *transa, const char *transb, const long int *m, 
 *                    const long int *n, const long int *k,              
 *                    const T *al, const T *a,                           
 *                    const long int *lda, const T *b,                   
 *                    const long int *ldb, const T *be, T *c,            
 *                    const long int *ldc)
 *
 * and call of style:
 *             gemm(transa,transb,m,n,k,al,a,lda,b,ldb,be,c,ldc)
 *
 * are all replaced by the unique definition call
 *
 *  template<class A, class B, class C, class Alpha, class Beta>
 *  void  mm(blas_status const& , A const& a, B const& a1, C& c,
 *           Alpha const& alpha = 1, Beta const& beta = 0)
 *
 * template<std::size_t TYPE = blas_types::general,
 *          char UPLO   = 'L',
 *          char TRANSA = 'N',
 *          char TRANSB = 'N',
 *          char DIAG   = 'N',
 *          char SIDE   = 'L'>
 * struct blas_status
 * {
 *   static const char uplo   = UPLO;
 *   static const char transa = TRANSA;
 *   static const char transb = TRANSB; 
 *   static const char diag   = DIAG;
 *   static const char side   = SIDE; 
 * };
 *
 * and call
 *   nt2::mm(nt2::blas_status<>(), a, b, r_);
 *
 * in the simpler case
 *
 * Note:
 * It must be noted that the fact that a matrix is triangular, symetric 
 * or hermitian can not be always known from the matrix shape. The routines do not
 * assume the matrix is really such but that only a part of the martix is significative
 * and the other never accessed.
 * For instance, a triangular matrix is never detected at compile time,
 * but it can be said to be such, which is different.
 **/
