/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_ALGEBRA_BLAS_MISC_LEVEL1_HPP_INCLUDED
#define NT2_TOOLBOX_ALGEBRA_BLAS_MISC_LEVEL1_HPP_INCLUDED

// SUBROUTINE _SWAP ( N,         X, INCX, Y, INCY )                              S, D, C, Z
// SUBROUTINE _COPY ( N,         X, INCX, Y, INCY )                              S, D, C, Z
// FUNCTION   _DOT  ( N,         X, INCX, Y, INCY )                              S, D, DS
// FUNCTION   _DOTU ( N,         X, INCX, Y, INCY )                              C, Z
// FUNCTION   _DOTC ( N,         X, INCX, Y, INCY )                              C, Z
// SUBROUTINE _AXPY ( N,  ALPHA, X, INCX, Y, INCY )                              S, D, C, Z
// FUNCTION   __DOT ( N,  ALPHA, X, INCX, Y, INCY )                              SDS
// SUBROUTINE _SCAL ( N,  ALPHA, X, INCX )                                       S, D, C, Z, CS, ZD
// FUNCTION   _NRM2 ( N,         X, INCX )                                       S, D, SC, DZ
// FUNCTION   _ASUM ( N,         X, INCX )                                       S, D, SC, DZ
// FUNCTION   I_AMAX( N,         X, INCX )                                       S, D, C, Z

namespace nt2
{
  namespace details
  {
    //swap
#define NT2_MISC(T, PREFIX, NAME)                  \
    inline void NAME(                              \
                      const long int *n,           \
                      T *x, const long int *incx,  \
                      T *y, const long int *incy,  \
                      )                            \
    {                                              \
      BOOST_PP_CAT(BOOST_PP_CAT(NAME,_)            \
                   (n, x, incx, y, incy);          \
    }                                              \

    NT2_MISC(double, d, swap)
    NT2_MISC(float,  s, swap)
    NT2_MISC(std::complex<double>, z, swap)
    NT2_MISC(std::complex<float>,  c, swap)
#undef NT2_MISC
      
    //copy 
#define NT2_MISC(T, PREFIX, NAME)                        \
      inline void NAME(                                  \
                     const long int *n,                  \
                      const T *x, const long int *incx,  \
                     T *y, const long int *incy,         \
                     )                                   \
      {                                                  \
        BOOST_PP_CAT(BOOST_PP_CAT(NAME,_)                \
             (n, x, incx, y, incy);                      \
      }                                                  \

    NT2_MISC(double, d, copy)
    NT2_MISC(float,  s, copy)
    NT2_MISC(std::complex<double>, z, copy)
    NT2_MISC(std::complex<float>,  c, copy)

    //copy dot dotc dotu
#define NT2_MISC(T, PREFIX, NAME, RET)             \
    inline RET NAME(                               \
                      const long int *n,           \
                      T *x, const long int *incx,  \
                      T *y, const long int *incy,  \
                      )                            \
    {                                              \
      BOOST_PP_CAT(BOOST_PP_CAT(NAME,_)            \
                   (n, x, incx, y, incy);          \
    }                                              \

      
    NT2_MISC(double, d, dot, double)
    NT2_MISC(float,  s, dot, float,)
    NT2_MISC(std::complex<double>, z, dotu, std::complex<double>)
    NT2_MISC(std::complex<float>,  c, dotu, std::complex<float>)
    NT2_MISC(std::complex<double>, z, dotc, std::complex<double>)
    NT2_MISC(std::complex<float>,  c, dotc, std::complex<float>)
      
#undef NT2_MISC

    //axpy
#define NT2_MISC(T, PREFIX, NAME)                  \
      inline T NAME(                               \
                      const long int *n,           \
                      const T *alpha,              \
                      T *x, const long int *incx,  \
                      )                            \
    {                                              \
      BOOST_PP_CAT(BOOST_PP_CAT(NAME,_)            \
                   (n, alpha, x, incx);            \
    }                                              \

    NT2_MISC(double, d, axpy)
    NT2_MISC(float,  s, axpy)
    NT2_MISC(std::complex<double>, z, axpy)
    NT2_MISC(std::complex<float>,  c, axpy)
      
#undef NT2_MISC

    //nrm2 asum i_amax
#define NT2_MISC(T, PREFIX, NAME, RET)             \
      inline RET NAME(                             \
                      const long int *n,           \
                      const T *alpha,              \
                      T *x, const long int *incx,  \
                      T *y, const long int *incy,  \
                      )                            \
    {                                              \
      BOOST_PP_CAT(BOOST_PP_CAT(NAME,_)            \
                   (n, alpha, x, incx, y, incy);   \
    }                                              \

    NT2_MISC(double, d, nrm2, double)
    NT2_MISC(float,  s, nrm2, float)
    NT2_MISC(std::complex<double>, z, nrm2, std::complex<double>)
    NT2_MISC(std::complex<float>,  c, nrm2, std::complex<float>)
    NT2_MISC(double, d, asum, double)
    NT2_MISC(float,  s, asum, float)
    NT2_MISC(std::complex<double>, z, asum, std::complex<double>)
    NT2_MISC(std::complex<float>,  c, asum, std::complex<float)
    NT2_MISC(double, d, i_amax, long int)
    NT2_MISC(float,  s, i_amax, long int)
    NT2_MISC(std::complex<double>, z, i_amax, long int)
    NT2_MISC(std::complex<float>,  c, i_amax, long int)
      
#undef NT2_MISC
      }
}


#endif

// /////////////////////////////////////////////////////////////////////////////
// End of misc_level1.hpp
// /////////////////////////////////////////////////////////////////////////////
