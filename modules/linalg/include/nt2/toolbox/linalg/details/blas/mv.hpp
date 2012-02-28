/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_DETAILS_BLAS_MV_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_DETAILS_BLAS_MV_HPP_INCLUDED
#include <nt2/toolbox/linalg/details/blas/blas2.hpp>

// all mv call: gemv symv hemv trmv gbmv sbmv hbmv tbmv etc.
// for data types float, double and related complex 

// _GEMV (        TRANS,      M, N,         ALPHA, A, LDA, X, INCX, BETA,  Y, INCY ) S, D, C, Z
// _GBMV (        TRANS,      M, N, KL, KU, ALPHA, A, LDA, X, INCX, BETA,  Y, INCY ) S, D, C, Z
// _HEMV ( UPLO,                 N,         ALPHA, A, LDA, X, INCX, BETA,  Y, INCY ) C, Z
// _HEMV ( UPLO,                 N,         ALPHA, A, LDA, X, INCX, BETA,  Y, INCY ) C, Z
// _SYMV ( UPLO,                 N,         ALPHA, A, LDA, X, INCX, BETA,  Y, INCY ) S, D
// _HBMV ( UPLO,                 N, K,      ALPHA, A, LDA, X, INCX, BETA,  Y, INCY ) C, Z
// _SBMV ( UPLO,                 N, K,      ALPHA, A, LDA, X, INCX, BETA,  Y, INCY ) S, D
// _HPMV ( UPLO,                 N,         ALPHA, AP,     X, INCX, BETA,  Y, INCY ) C, Z
// _SPMV ( UPLO,                 N,         ALPHA, AP,     X, INCX, BETA,  Y, INCY ) S, D
// _TPMV ( UPLO, TRANS, DIAG,    N,                AP,     X, INCX )                 S, D, C, Z
// _TRMV ( UPLO, TRANS, DIAG,    N,                A, LDA, X, INCX )                 S, D, C, Z
// _TBMV ( UPLO, TRANS, DIAG,    N, K,             A, LDA, X, INCX )                 S, D, C, Z

namespace nt2
{
  namespace details
  {
#define NT2_MV(T, PREFIX)                                               \
    inline void gemv(const char *ta, const long int *m,                 \
                     const long int *n,                                 \
                     const T *al,                                       \
                     const T *a, const long int *lda,                   \
                     const T *x, const long int *incx,                  \
                     const T *be,                                       \
                     T *y, const long int *incy)                        \
    {                                                                   \
      BOOST_PP_CAT(PREFIX,BOOST_PP_CAT(gemv,_))                         \
        (ta,m,n,al,a,lda,x,incx,be,y,incy);                             \
    }                                                                   \
        
    NT2_MV(double, d)
    NT2_MV(float,  s)
    NT2_MV(std::complex<double>, z)
    NT2_MV(std::complex<float>, c)

#undef NT2_MV

#define NT2_MV(T, PREFIX)                                               \
    inline void gbmv(const char *ta, const long int *m,                 \
                     const long int *n,                                 \
                     const long int *kl, const long int *ku,            \
                     const T *al,                                       \
                     const T *a, const long int *lda,                   \
                     const T *x, const long int *incx,                  \
                     const T *be,                                       \
                     T *y, const long int *incy)                        \
    {                                                                   \
      BOOST_PP_CAT(PREFIX,BOOST_PP_CAT(gbmv,_))                         \
        (ta,m,n,kl,ku,al,a,lda,x,incx,be,y,incy);                       \
    }                                                                   \
        
    NT2_MV(double, d)
    NT2_MV(float,  s)
    NT2_MV(std::complex<double>, z)
    NT2_MV(std::complex<float>, c)

#undef NT2_MV
      
#define NT2_MV(T, PREFIX, PREFIX2)                                      \
      inline void  BOOST_PP_CAT(PREFIX2,mv)(const char *uplo,           \
                     const long int *n,                                 \
                     const T *al,                                       \
                     const T *a, const long int *lda,                   \
                     const T *x, const long int *incx,                  \
                     const T *be,                                       \
                     T *y, const long int *incy)                        \
    {                                                                   \
      BOOST_PP_CAT(PREFIX,BOOST_PP_CAT(BOOST_PP_CAT(PREFIX2,mv),_))     \
        (ta,n,al,a,lda,x,incx,be,y,incy);                               \
    }                                                                   \
        
    NT2_MV(std::complex<double>, z, he)
    NT2_MV(std::complex<float>, c, he)
    NT2_MV(std::complex<double>, z, sy)
    NT2_MV(std::complex<float>, c, sy)
    NT2_MV(double, d, sy)
    NT2_MV(float,  s, sy)

#undef NT2_MV

#define NT2_MV(T, PREFIX, PREFIX2)                                      \
      inline void BOOST_PP_CAT(PREFIX2,mv)(const char *uplo,            \
                       const long int *n, const long int *k,            \
                       const T *al,                                     \
                       const T *a, const long int *lda,                 \
                       const T *x, const long int *incx,                \
                       const T *be,                                     \
                       T *y, const long int *incy)                      \
    {                                                                   \
      BOOST_PP_CAT(PREFIX,BOOST_PP_CAT(BOOST_PP_CAT(PREFIX2,mv),_))     \
        (ta,n,k,al,a,lda,x,incx,be,y,incy);                             \
    }                                                                   \
        
    NT2_MV(std::complex<double>, z, hb)
    NT2_MV(std::complex<float>, c, hb)
    NT2_MV(std::complex<double>, z, sb)
    NT2_MV(std::complex<float>, c, sb)
    NT2_MV(double, d, sb)
    NT2_MV(float,  s, sb)

#undef NT2_MV

#define NT2_MV(T, PREFIX)                                               \
      inline void BOOST_PP_CAT(PREFIX2,mv)(const char *uplo,            \
                     const long int *n,                                 \
                     const T *al,                                       \
                     const T *ap,                                       \
                     const T *x, const long int *incx,                  \
                     const T *be,                                       \
                     T *y, const long int *incy)                        \
    {                                                                   \
      BOOST_PP_CAT(PREFIX,BOOST_PP_CAT(hpmv,_))                         \
        (ta,n,al,ap,x,incx,be,y,incy);                                  \
    }                                                                   \
        
    NT2_MV(std::complex<double>, z, hp)
    NT2_MV(std::complex<float>, c, hp)
    NT2_MV(std::complex<double>, z, sp)
    NT2_MV(std::complex<float>, c, sp)
    NT2_MV(double, d, sp)
    NT2_MV(float,  s, sp)

#undef NT2_MV

#define NT2_MV(T, PREFIX)                                               \
    inline void tpmv(const char *uplo, const char *trans,               \
                     const char *diag,                                  \
                     const long int *n,                                 \
                     const T *ap,                                       \
                     T *x, const long int *incx)                        \
    {                                                                   \
      BOOST_PP_CAT(PREFIX,BOOST_PP_CAT(tpmv,_))                         \
        (uplo,trans,diag,n,ap,x,incx);                               \
    }                                                                   \
        
      NT2_MV(double, d)
      NT2_MV(float,  s)
      NT2_MV(std::complex<double>, z)
      NT2_MV(std::complex<float>, c)
      
#undef NT2_MV

#define NT2_MV(T, PREFIX)                                               \
    inline void trmv(const char *uplo, const char *trans,               \
                     const char *diag,                                  \
                     const long int *n,                                 \
                     const T *a, const long int *lda                    \
                     T *x, const long int *incx)                        \
    {                                                                   \
      BOOST_PP_CAT(PREFIX,BOOST_PP_CAT(tpmv,_))                         \
        (uplo,trans,diag,n,a,lda,x,incx);                               \
    }                                                                   \
        
      NT2_MV(double, d)
      NT2_MV(float,  s)
      NT2_MV(std::complex<double>, z)
      NT2_MV(std::complex<float>, c)
      
#undef NT2_MV
           
      
#define NT2_MV(T, PREFIX)                                               \
    inline void tbmv(const char *uplo, const char *trans,               \
                     const char *diag,                                  \
                     const long int *n,  const long int *k,             \
                     const T *a,                                        \
                     const long int *lda,                               \
                     T *x, const long int *incx)                        \
    {                                                                   \
      BOOST_PP_CAT(PREFIX,BOOST_PP_CAT(trmv,_))                         \
        (uplo,trans,diag,n,k, a,lda,x,incx);                            \
    }                                                                   \
        
      NT2_MV(double, d)
      NT2_MV(float,  s)
      NT2_MV(std::complex<double>, z)
      NT2_MV(std::complex<float>, c)
      
#undef NT2_MV

        
  }  
}


#endif

// /////////////////////////////////////////////////////////////////////////////
// End of mv.hpp
// /////////////////////////////////////////////////////////////////////////////
