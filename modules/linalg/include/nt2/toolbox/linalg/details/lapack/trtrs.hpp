/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_DETAILS_LAPACK_TRTRS_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_DETAILS_LAPACK_TRTRS_HPP_INCLUDED
// trtrs

namespace nt2
{
  namespace details
  {
    extern "C"
    {
#define NT2_COMPLEX void
      void NT2_F77NAME(ctrtrs)(const char* uplo, const char* trans, const char* diag,
                               const long int* n, const long int* nrhs, const NT2_COMPLEX* a, const long int* lda,
                               NT2_COMPLEX* b, const long int* ldb, long int* info);
      void NT2_F77NAME(dtrtrs)(const char* uplo, const char* trans, const char* diag,
                               const long int* n, const long int* nrhs, const double* a, const long int* lda,
                               double* b, const long int* ldb, long int* info);
      void NT2_F77NAME(strtrs)(const char* uplo, const char* trans, const char* diag,
                               const long int* n, const long int* nrhs, const float* a, const long int* lda,
                               float* b, const long int* ldb, long int* info);
      void NT2_F77NAME(ztrtrs)(const char* uplo, const char* trans, const char* diag,
                               const long int* n, const long int* nrhs, const NT2_COMPLEX* a, const long int* lda,
                               NT2_COMPLEX* b, const long int* ldb, long int* info);
#undef NT2_COMPLEX
    }
    
#define NT2_TRTRS(NAME, T)                      \
    inline void trtrs(const char* uplo,         \
                      const char* trans,        \
                      const char* diag,         \
                      const long int* n,        \
                      const long int* nrhs,     \
                      const T* a,               \
                      const long int* lda,      \
                      T* b,                     \
                      const long int* ldb,      \
                      long int* info)           \
    {                                           \
      NT2_F77NAME( NAME )(uplo,trans,diag,      \
                          n,nrhs,a,             \
                          lda, b,ldb,info);     \
    }                                           \


    NT2_TRTRS(strtrs, float)
    NT2_TRTRS(dtrtrs, double)
    NT2_TRTRS(ctrtrs, std::complex<float>)
    NT2_TRTRS(ztrtrs, std::complex<double>)
      
#undef NT2_TRTRS
  }
}

#endif

// /////////////////////////////////////////////////////////////////////////////
// End of trs.hpp
// /////////////////////////////////////////////////////////////////////////////
