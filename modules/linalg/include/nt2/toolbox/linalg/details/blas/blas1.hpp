//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II       
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI       
//                                                                              
//          Distributed under the Boost Software License, Version 1.0.          
//                 See accompanying file LICENSE.txt or copy at                 
//                     http://www.boost.org/LICENSE_1_0.txt                     
//==============================================================================
#ifndef NT2_TOOLBOX_LINALG_DETAILS_BLAS_BLAS1_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_DETAILS_BLAS_BLAS1_HPP_INCLUDED

#include <nt2/toolbox/linalg/details/utility/f77_wrapper.hpp>

extern "C"
{
  // Real, double precision  
  double NT2_F77NAME(dasum)(const long int *n, const double *dx, const long int *incx);
    
  void NT2_F77NAME(daxpy)(const long int *n, const double *da, const double *dx, 
                      const long int *incx, double *dy, const long int *incy);
  
  void NT2_F77NAME(dcopy)(const long int *n, double *dx, const long int *incx, 
                      double *dy, const long int *incy);
  
  
  double NT2_F77NAME(ddot)(const long int *n, const double *dx, const long int *incx, 
                       const double *dy, const long int *incy);
  
  double NT2_F77NAME(dnrm2)(const long int *n, const double *dx, const long int *incx); 
  
  void NT2_F77NAME(drot)(const long int *n, double *dx, const long int *incx, 
                     double *dy, const long int *incy, const double *c, 
                     const double *s);
  
  void NT2_F77NAME(drotg)(double *da, double *db, double *c, double *s);
  
  void NT2_F77NAME(dscal)(const long int *n, double *da, double *dx, 
                      const long int *incx);
  
  void NT2_F77NAME(dswap)(const long int *n, double *dx, const long int *incx, 
                      double *dy, const long int *incy);
  
  long int NT2_F77NAME(idamax)(const long int *n, const double *dx, 
                           const long int *incx);
  
  // Real, single precision
  
  float NT2_F77NAME(sasum)(const long int *n, const float *dx, const long int *incx);
  
  
  void NT2_F77NAME(saxpy)(const long int *n, const float *da, const float *dx, 
                      const long int *incx, float *dy, const long int *incy);
  
  void NT2_F77NAME(scopy)(const long int *n, float *dx, const long int *incx, float *dy, 
                      const long int *incy);
  
  
  float NT2_F77NAME(sdot)(const long int *n, const float *dx, const long int *incx, 
                      const float *dy, const long int *incy);
  
  float NT2_F77NAME(snrm2)(const long int *n, const float *dx, const long int *incx); 
  
  void NT2_F77NAME(srot)(const long int *n, float *dx, const long int *incx, float *dy, 
                     const long int *incy, const float *c, const float *s);
  
  void NT2_F77NAME(srotg)(float *da, float *db, float *c, float *s);
  
  void NT2_F77NAME(sscal)(const long int *n, float *da, float *dx, const long int *incx);
  
  void NT2_F77NAME(sswap)(const long int *n, float *dx, const long int *incx, float *dy, 
                      const long int *incy);
  
  long int NT2_F77NAME(isamax)(const long int *n, const float *dx, const long int *incx);
  
#define NT2_WRAP_COMPLEX void
  
  
  //complex < float >
  
  float NT2_F77NAME(cdotc)(NT2_WRAP_COMPLEX *c, const long int *n, 
                       const NT2_WRAP_COMPLEX *cx, const long int *incx, 
                       const NT2_WRAP_COMPLEX *cy, const long int *incy);
  
  float NT2_F77NAME(cdotu)(NT2_WRAP_COMPLEX *c, const long int *n, 
                       const NT2_WRAP_COMPLEX *cx, const long int *incx, 
                       const NT2_WRAP_COMPLEX *cy, const long int *incy);
  
  void NT2_F77NAME(caxpy)(const long int *n, const NT2_WRAP_COMPLEX *da, 
                      const NT2_WRAP_COMPLEX *dx, 
                      const long int *incx, NT2_WRAP_COMPLEX *dy, 
                      const long int *incy);
  
  void NT2_F77NAME(ccopy)(const long int *n, NT2_WRAP_COMPLEX *dx, const long int *incx, 
                      NT2_WRAP_COMPLEX *dy, const long int *incy);
  
  float  NT2_F77NAME(scasum)(const long int *n, const NT2_WRAP_COMPLEX *dx, 
                         const long int *incx);
  
  float  NT2_F77NAME(scnrm2)(const long int *n, const NT2_WRAP_COMPLEX *dx, 
                         const long int *incx); 
  
  void NT2_F77NAME(cdscal)(const long int *n, const float *da, NT2_WRAP_COMPLEX *dx, 
                       const long int *incx);
  
  void NT2_F77NAME(cscal)(const long int *n, const NT2_WRAP_COMPLEX *da, 
                      NT2_WRAP_COMPLEX *dx, const long int *incx);
  
  long int NT2_F77NAME(icamax)(const long int *n, const NT2_WRAP_COMPLEX *dx, 
                           const long int *incx);
  
  void NT2_F77NAME(cswap)(const long int *n, NT2_WRAP_COMPLEX *dx, const long int *incx, 
                      NT2_WRAP_COMPLEX *dy, long int *incy);
  
  //complex < double > 
  
  double NT2_F77NAME(zdotc)(NT2_WRAP_COMPLEX *c, const long int *n, 
                        const NT2_WRAP_COMPLEX *cx, 
                        const long int *incx, const NT2_WRAP_COMPLEX *cy, 
                        const long int *incy);
  
  double NT2_F77NAME(zdotu)(NT2_WRAP_COMPLEX *c, const long int *n, 
                        const NT2_WRAP_COMPLEX *cx, const long int *incx, 
                        const NT2_WRAP_COMPLEX *cy, const long int *incy);
  
  void NT2_F77NAME(zaxpy)(const long int *n, const NT2_WRAP_COMPLEX *da, 
                      const NT2_WRAP_COMPLEX *dx, 
                      const long int *incx, NT2_WRAP_COMPLEX *dy, 
                      const long int *incy);
  
  void NT2_F77NAME(zcopy)(const long int *n, NT2_WRAP_COMPLEX *dx, const long int *incx, 
                      NT2_WRAP_COMPLEX *dy, const long int *incy);
  
  double  NT2_F77NAME(dzasum)(const long int *n, const NT2_WRAP_COMPLEX *dx, 
                          const long int *incx);
  
  double  NT2_F77NAME(dznrm2)(const long int *n, const NT2_WRAP_COMPLEX *dx, 
                          const long int *incx); 
  
  void NT2_F77NAME(zdscal)(const long int *n, const double *da, NT2_WRAP_COMPLEX *dx, 
                       const long int *incx);
  
  void NT2_F77NAME(zscal)(const long int *n, const NT2_WRAP_COMPLEX *da, 
                      NT2_WRAP_COMPLEX *dx, const long int *incx);
  
  long int NT2_F77NAME(izamax)(const long int *n, const NT2_WRAP_COMPLEX *dx, 
                           const long int *incx);
  
  void NT2_F77NAME(zswap)(const long int *n, NT2_WRAP_COMPLEX *dx, const long int *incx, 
                      NT2_WRAP_COMPLEX *dy, long int *incy);
  
#undef NT2_WRAP_COMPLEX
  
}

#endif
