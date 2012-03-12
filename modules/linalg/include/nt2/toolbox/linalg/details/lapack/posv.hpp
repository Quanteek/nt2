/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2009-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_DETAILS_LAPACK_POSV_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_DETAILS_LAPACK_POSV_HPP_INCLUDED
#include <nt2/toolbox/linalg/details/utility/f77_wrapper.hpp>
#include <nt2/toolbox/linalg/details/lapack/workspace.hpp>

namespace nt2
{
  namespace details
  {
    extern "C"
    {
#define NT2_COMPLEX void
      void NT2_F77NAME(cposv)(const char* uplo, const long int* n, const long int* nrhs,
                              COMPLEX* a, const long int* lda, COMPLEX* b, const long int* ldb,
                              long int* info);
      void NT2_F77NAME(dposv)(const char* uplo, const long int* n, const long int* nrhs,
                              double* a, const long int* lda, double* b, const long int* ldb,
                              long int* info);
      void NT2_F77NAME(sposv)(const char* uplo, const long int* n, const long int* nrhs,
                              float* a, const long int* lda, float* b, const long int* ldb,
                              long int* info);
      void NT2_F77NAME(zposv)(const char* uplo, const long int* n, const long int* nrhs,
                              COMPLEX* a, const long int* lda, COMPLEX* b, const long int* ldb,
                              long int* info);
#undef NT2_COMPLEX
    }

#define NT2_POSV(NAME, T)                       \
    inline void poposv(const char* uplo,        \
                       const long int* n,       \
                       const long int* nrhs,    \
                       T* a,                    \
                       const long int* lda,     \
                       T* b,                    \
                       const long int* ldb,     \
                       long int* info)          \
      {                                         \
        NT2_F77NAME( NAME )(uplo,n,nrhs,        \
                            a,lda,b,ldb,info);  \
      }                                         \
          
    NT2_POPOSV(sposv, float)
    NT2_POPOSV(dposv, double)
    NT2_POPOSV(cposv, std::complex<float>)
    NT2_POPOSV(zposv, std::complex<double>)

#undef NT2_POSV


  }
}


#endif

