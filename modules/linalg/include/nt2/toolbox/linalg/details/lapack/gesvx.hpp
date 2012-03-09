/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_DETAILS_LAPACK_GESVX_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_DETAILS_LAPACK_GESVX_HPP_INCLUDED
namespace nt2
{
  namespace details
  {
    extern "C"
    {
#define NT2_COMPLEX void
      void NT2_F77NAME(cgesvx)(const char* fact, const char* trans, const long int* n, const long int* nrhs,
                               const NT2_COMPLEX* a, const long int* lda, NT2_COMPLEX* af, const long int* ldaf, long int* ipiv,
                               char* equed, const float* r, const float* c, NT2_COMPLEX* b, const long int* ldb,
                               NT2_COMPLEX* x, const long int* ldx, float* rcond, float* ferr, float* berr, NT2_COMPLEX* work,
                               float* rwork, long int* info);
      void NT2_F77NAME(dgesvx)(const char* fact, const char* trans, const long int* n, const long int* nrhs,
                               const double* a, const long int* lda, double* af, const long int* ldaf, long int* ipiv,
                               char* equed, const double* r, const double* c, double* b, const long int* ldb,
                               double* x, const long int* ldx, double* rcond, double* ferr, double* berr, double* work,
                               long int* iwork, long int* info);
      void NT2_F77NAME(sgesvx)(const char* fact, const char* trans, const long int* n, const long int* nrhs,
                               const float* a, const long int* lda, float* af, const long int* ldaf, long int* ipiv,
                               char* equed, const float* r, const float* c, float* b, const long int* ldb,
                               float* x, const long int* ldx, float* rcond, float* ferr, float* berr, float* work,
                               long int* iwork, long int* info);
      void NT2_F77NAME(zgesvx)(const char* fact, const char* trans, const long int* n, const long int* nrhs,
                               const NT2_COMPLEX* a, const long int* lda, NT2_COMPLEX* af, const long int* ldaf, long int* ipiv,
                               char* equed, const double* r, const double* c, NT2_COMPLEX* b, const long int* ldb,
                               NT2_COMPLEX* x, const long int* ldx, double* rcond, double* ferr, double* berr, NT2_COMPLEX* work,
                               double* rwork, long int* info);
#undef NT2_COMPLEX
    }
    
#define NT2_GESVX(NAME, T, TBASE)                                       \
    inline void gesvx(const char* fact,                                 \
                      const char* trans,                                \
                      const long int* n,                                \
                      const long int* nrhs,                             \
                      const T* a,                                       \
                      const long int* lda,                              \
                      T* af,                                            \
                      const long int* ldaf,                             \
                      long int* ipiv,                                   \
                      char* equed,                                      \
                      const TBASE* r,                                   \
                      const TBASE* c,                                   \
                      T* b,                                             \
                      const long int* ldb,                              \
                      T* x,                                             \
                      const long int* ldx,                              \
                      TBASE* rcond,                                     \
                      TBASE* ferr,                                      \
                      TBASE* berr,                                      \
                      long int* info,                                   \
                      workspace<T> & w)                                 \
    {                                                                   \
      w.resizerw(2**n);                                                 \
      w.resizew(2**n);                                                  \
      NT2_F77NAME( NAME )(fact, trans, n, nrhs, a, lda, af, ldaf,       \
                          ipiv, equed, r, c, b, ldb, x, ldx,            \
                          rcond, ferr, berr, w.getw(),                  \
                          w.getrw(), info);                             \
    }                                                                   \
    inline void gesvx(const char* fact,                             \
                      const char* trans,                                \
                      const long int* n,                                \
                      const long int* nrhs,                             \
                      const T* a,                                       \
                      const long int* lda,                              \
                      T* af,                                            \
                      const long int* ldaf,                             \
                      long int* ipiv,                                   \
                      char* equed,                                      \
                      const TBASE* r,                                   \
                      const TBASE* c,                                   \
                      T* b,                                             \
                      const long int* ldb,                              \
                      T* x,                                             \
                      const long int* ldx,                              \
                      TBASE* rcond,                                     \
                      TBASE* ferr,                                      \
                      TBASE* berr,                                      \
                      long int* info)                                   \
    {                                                                   \
      workspace<T> w;                                                   \
      gesvx(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv,               \
            equed, r, c, b, ldb, x, ldx,                                \
            rcond, ferr, berr, info, w);                                \
    }                                                                   \
        
    NT2_GESVX(cgesvx, float,  float)
    NT2_GESVX(zgesvx, double, double)
    NT2_GESVX(cgesvx, std::complex<float>,  float)
    NT2_GESVX(zgesvx, std::complex<double>, double)
      
#undef NT2_GESVX
  }
}

#endif

// /////////////////////////////////////////////////////////////////////////////
// End of svx.hpp
// /////////////////////////////////////////////////////////////////////////////
