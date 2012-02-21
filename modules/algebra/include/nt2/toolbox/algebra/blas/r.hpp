/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_ALGEBRA_BLAS_R_HPP_INCLUDED
#define NT2_TOOLBOX_ALGEBRA_BLAS_R_HPP_INCLUDED

// _GER  (                    M, N, ALPHA, X, INCX, Y, INCY, A, LDA ) S, D
// _GERU (                    M, N, ALPHA, X, INCX, Y, INCY, A, LDA ) C, Z (TRANSPOSE)
// _GERC (                    M, N, ALPHA, X, INCX, Y, INCY, A, LDA ) C, Z (TRANSCONJUGATE)
// _HER  ( UPLO,                 N, ALPHA, X, INCX,          A, LDA ) C, Z
// _SYR  ( UPLO,                 N, ALPHA, X, INCX,          A, LDA ) S, D
// _HER2 ( UPLO,                 N, ALPHA, X, INCX, Y, INCY, A, LDA ) C, Z
// _SYR2 ( UPLO,                 N, ALPHA, X, INCX, Y, INCY, A, LDA ) S, D
// _SYR2 ( UPLO,                 N, ALPHA, X, INCX, Y, INCY, A, LDA ) S, D
// _HPR2 ( UPLO,                 N, ALPHA, X, INCX, Y, INCY, AP )     C, Z
// _SPR2 ( UPLO,                 N, ALPHA, X, INCX, Y, INCY, AP )     S, D
// _SPR  ( UPLO,                 N, ALPHA, X, INCX,          AP )     S, D
// _HPR  ( UPLO,                 N, ALPHA, X, INCX,          AP )     C, Z

// gerx A <- {t,c}x*y+A      with x =  u, c or nothing in which case x == y
// herx A <- hx*y+hy*x+A     with x = 2 or nothing in which case x == y
// syrx A <- tx*y+ty*x+A     with x = 2 or nothing in which case x == y

namespace nt2
{
  namespace details
  {
    //ger geru gerc
#define NT2_R(T, PREFIX, SUFFIX)                         \
    inline void BOOST_PP_CAT(ger,SUFFIX)(                \
                     const long int *m,                  \
                     const long int *n,                  \
                     const T *al,                        \
                     const T *x, const long int *incx,   \
                     const T *y, const long int *incy    \
                     const T *a, const long int *lda     \
                       )                                 \
    {                                                    \
      BOOST_PP_CAT(BOOST_PP_CAT(ger,SUFFIX),_))          \
        (ta,m,n,al,x,incx,y,incy,a,lda);                 \
    }                                                    \
        
    NT2_R(double, d, )
    NT2_R(float,  s, )
    NT2_R(std::complex<double>, z, u)
    NT2_R(std::complex<float>,  c, u)
    NT2_R(std::complex<double>, z, c)
    NT2_R(std::complex<float>,  c, c)

#undef NT2_R

    //her syr
#define NT2_R(T, PREFIX, PREFIX2)                                   \
      inline void BOOST_PP_CAT(PREFIX2,r)(                          \
                    const char * uplo,                              \
                    const long int *m,                              \
                    const long int *n,                              \
                    const T *al,                                    \
                    const T *x, const long int *incx,               \
                    const T *be,                                    \
                    T *a, const long int *lda )                     \
    {                                                               \
      BOOST_PP_CAT(BOOST_PP_CAT(PREFIX,BOOST_PP_CAT(PREFIX2,r),_))  \
        (ta,m,n,al,x,incx,y,incy,a,lda);                            \
    }                                                               \
        
    NT2_R(std::complex<double>, z, he)
    NT2_R(std::complex<float>,  c, he)
    NT2_R(std::complex<double>, s, sy)
    NT2_R(std::complex<float>,  d, sy)
#undef NT2_R
      
    //hpr spr
#define NT2_R(T, PREFIX, PREFIX2)                                   \
      inline void BOOST_PP_CAT(PREFIX2,r)(                          \
                    const char * uplo,                              \
                    const long int *m,                              \
                    const long int *n,                              \
                    const T *al,                                    \
                    const T *x, const long int *incx,               \
                    const T *be,                                    \
                    T *ap )                                         \
    {                                                               \
      BOOST_PP_CAT(BOOST_PP_CAT(PREFIX,BOOST_PP_CAT(PREFIX2,r),_))  \
        (ta,m,n,al,x,incx,y,incy,ap);                               \
    }                                                               \
        
    NT2_R(std::complex<double>, z, hp)
    NT2_R(std::complex<float>,  c, hp)
    NT2_R(std::complex<double>, s, sp)
    NT2_R(std::complex<float>,  d, sp)
#undef NT2_R
  }

//her2 syr2
#define NT2_R(T, PREFIX, PREFIX2)                                   \
  inline void BOOST_PP_CAT(PREFIX2,r2)(                             \
                    const char * uplo,                              \
                    const long int *m,                              \
                    const long int *n,                              \
                    const T *al,                                    \
                    const T *x, const long int *incx,               \
                    const T *y, const long int *incy,               \
                    const T *be,                                    \
                    T *a, const long int *lda)                      \
    {                                                               \
      BOOST_PP_CAT(BOOST_PP_CAT(PREFIX,BOOST_PP_CAT(PREFIX2,r2),_)) \
        (ta,m,n,al,x,incx,y,incy,a, lda);                           \
    }                                                               \
        
    NT2_R(std::complex<double>, z, he)
    NT2_R(std::complex<float>,  c, he)
    NT2_R(std::complex<double>, s, sy)
    NT2_R(std::complex<float>,  d, sy)
#undef NT2_R

    //hpr spr
#define NT2_R(T, PREFIX, PREFIX2)                                   \
      inline void BOOST_PP_CAT(PREFIX2,r)(                          \
                    const char * uplo,                              \
                    const long int *m,                              \
                    const long int *n,                              \
                    const T *al,                                    \
                    const T *x, const long int *incx,               \
                    const T *be,                                    \
                    T *ap )                                         \
    {                                                               \
      BOOST_PP_CAT(BOOST_PP_CAT(PREFIX,BOOST_PP_CAT(PREFIX2,r2),_)) \
        (ta,m,n,al,x,incx,ap);                                      \
    }                                                               \
        
    NT2_R(std::complex<double>, z, hp)
    NT2_R(std::complex<float>,  c, hp)
    NT2_R(std::complex<double>, s, sp)
    NT2_R(std::complex<float>,  d, sp)
#undef NT2_R

      
  }
}


#endif

// /////////////////////////////////////////////////////////////////////////////
// End of r.hpp
// /////////////////////////////////////////////////////////////////////////////
