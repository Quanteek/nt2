#error "This file is for documentation purpose only."
/*!
 * \ingroup toolbox
 * \defgroup lapack lapack functions
 * \brief Defines lapack functions
 **/
/*!
 * \ingroup lapack
 * \defgroup lapack
 * \brief Defines lapack
 *
 * The aim of this part is to encapsulate the lapack functions
 * in a generic manner suited to C++
 *
 * rationale
 * lapack functions consits of name and parameters
 *
 * Lapack names:
 * lapack function have names consisting in 3 parts:
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
 * Lapack parameters
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
 * of performed operation prefixed with b_ (as lapack) and will use much less parameters than the original ones
 *
 * * The data pointers will be nt2 containers which are aware of their shape and of their
 *   dimensions and extents,  suppressing the need of passing dimensions as supplementary long ints
 *   parameters as well as storage information in the function name.
 *
 * * The char properties can be subsumed in a struct containing only static constants allowing to
 *   know the properties at compile time (but disallowing to chose to change dynamically the property, 
 *   we think this is rarely needed and can be replaced by testing and branching if not avoidable)
 *
 * ---------------------------------------------------------------------------------------------------
 *
 * So for instance, the functions gesv,  gbsv, hesv, hbsv, sysv, sbsv etc.
 * which all solve in various cases
 *                               a*b = c
 * with 40 FORTRAN distinct definition calls analog to the following C call
 * for each type and shape
 *    void gesv(const char *transa,
 *                    const long int *n, const long int *nrsh,              
 *                    T *a, const long int *lda,
 *                    T *b, const long int *ldb,
 *                    T *c, const long int *ldc,
 *                    const long int *info)
 *
 * and call of style:
 *             gesv(transa,n,nrsh,a,lda,b,ldb,c,ldc,info)
 *
 * are all replaced by the unique definition call (and few variants for returning more info)
 *
 *  template<class status, class A, class B, class C>
 *  bool sv(status const& , A const& a, B const& a1, C& c)
 *
 * template<std::size_t TYPE = lapack_types::general,
 *          char UPLO   = 'L',
 *          char TRANSA = 'N',
 *          char TRANSB = 'N',
 *          char DIAG   = 'N',
 *          char SIDE   = 'L'>
 * struct lapack_status
 * {
 *   static const char uplo   = UPLO;
 *   static const char transa = TRANSA;
 *   static const char transb = TRANSB; 
 *   static const char diag   = DIAG;
 *   static const char side   = SIDE; 
 * };
 *
 * and call
 *   nt2::sv(St(), a, b, r_);
 *
 * Where St is the corresponding status struct type
 *
 * Notes:
 * It must be noted that the fact that a matrix is triangular, symetric 
 * or hermitian can not be always known from the matrix shape. The routines do not
 * assume the matrix is really such but that only a part of the martix is significative
 * and the other never accessed.
 * For instance, a triangular matrix is never detected at compile time,
 * but it can be said to be such, which is different.
 
 * These are the type/form of matrices supported by lapack in alphabetical order
 * those prefixed with # where considered in nt2 v2
 * #BD  bidiagonal
 * #DI  diagonal
 * #GB  general band
 * GE   general (i.e., unsymmetric, in some cases rectangular)
 * GG   general matrices, generalized problem (i.e., a pair of general matrices)
 * #GT  general tridiagonal
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
 * #TR  triangular (or in some cases quasi-triangular)
 * TZ   trapezoidal
 * #UN  (complex) unitary
 * UP   (complex) unitary, packed storage

 Must be added types to hold the diagonal,  bidiagonal and tridiagonal matrices

 In Lapack (ecerpt from the lapack doc
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

 TO DO
 J'en viens à la conclusion qu'une matrice devrait savoir si elle est stockée normalement ou comme transposée
 typiquement gemm sait calculer A°*B° ou le ° est N T ou C
 mais symm ou A est symtrique ne sait calculer que AB et BA
 mais comme BA est le transposé de AB'
 symm sait calculer AB' mais sous forme transposée...
 en général on n'a pas besoin d'avoir une matrice dans un sens ou un autre
 mais juste de savoir qu'elle est dans un sens ou un autre et en profiter pour faire le bon
 calcul: on n'aurait jamais à vraiment transposer de matrices sauf si on veut sauf si on veut
 forcer la chose pour passer la matrices à des routines externes moins intelligentes.

 Cependant reste un pb c'est le dimensionnement des matrices de retour si on fait C = A*B'
 avec A mxm et B nxn 
 1 ) si on détecte que A est symetrique,  il faudrait calculer dans C la matrice BA qui est
 la transposé du résultat et dire qu'elle est stockée en mode transpose et ses dimensions seraient
 nxm
 2 )si on détecte que A est generale, il faudrait calculer dans C la matrice AB' effective qui est
 de dimension mxn
 Ca pose un pb + gros encore si c'est C =  aAB'+bC qu'on veut faire
 parce que là on ne veut pas transposer C !
 peut faut-il decider de transposer ou pas B suivant les dimensions de C ?! 

 On pourrait imaginer avoir un truc statique qui dise general_ ou general_do_not_specialize_
 le general_do_not_specialize_ se passerait des tests qui permettrait au run time
 d'utiliser des methodes statiquement definies
 pour des formes comme symetric. Il est impossible de savoir sur la forme d'une matrice generale,
 quelle est la nature de la matrice et il faut donc rajouter explicitement
 un parametre (statique) de nature (symetric etc.)

 je pense que le nom du  tag rectangular_ pour les tables nt2 est mal choisi ce
 devrait etre general_  (et  general_do_not_specialize_)
 ce qui serait bien plus compatible avec lapack et aurait bien plus de sens
 puisque de toute facon les tables sont rectangulaires et les schémas de stockage
 band et packed le sont aussi

 il y a 2 notions
   l'ordre dans le stockage general band packed
   la forme  symetric triangulaire etc
   
 mais la donnée d'une matrice triangulaire peut représenter une matrice symetrique si le stockage est general
 il y a une difference entre etre triangulaire  et etre considérée comme triangulaire
 
 le packed storage des matrices triangulaires et symetriques peut etre identique
 il faut qu'a un niveau quelconque les 2 notions soient distinguées.
 Pour un meme storage general a(i,j) peut designer 7 choses differentes
 suivant que la matrice est quelconque symetrique hermitienne ou triangulaire
 et que la routine est appelée avec UPLO=u ou l

 je voudrais qu'on accorde les violons
   soit table container etc. st capables de fournir tout ça
   soit à quoi faut-il se limiter
 
**/
