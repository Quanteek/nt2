/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_ALGEBRA_LAPACK_L_SV_C_HPP_INCLUDED
#define NT2_TOOLBOX_ALGEBRA_LAPACK_L_SV_C_HPP_INCLUDED
#include <nt2/toolbox/algebra/blas/f77_wrapper.hpp>

extern "C"
{
#define NT2_WRAP_COMPLEX void
  
 void NT2_F77NAME(cgbsv)(const long int* n, const long int* kl, const long int* ku, const long int* nrhs, NT2_WRAP_COMPLEX* ab, const long int* ldab, long int* ipiv, NT2_WRAP_COMPLEX* b, const long int* ldb, long int* info);
 void NT2_F77NAME(cgesv)(const long int* n, const long int* nrhs, NT2_WRAP_COMPLEX* a, const long int* lda, long int* ipiv, NT2_WRAP_COMPLEX* b, const long int* ldb, long int* info);
 void NT2_F77NAME(cgtsv)(const long int* n, const long int* nrhs, NT2_WRAP_COMPLEX* dl, NT2_WRAP_COMPLEX* d, NT2_WRAP_COMPLEX* du, NT2_WRAP_COMPLEX* b, const long int* ldb, long int* info);
 void NT2_F77NAME(chesv)(const char* uplo, const long int* n, const long int* nrhs, NT2_WRAP_COMPLEX* a, const long int* lda, long int* ipiv, NT2_WRAP_COMPLEX* b, const long int* ldb, NT2_WRAP_COMPLEX* work, const long int* lwork, long int* info);
 void NT2_F77NAME(chpsv)(const char* uplo, const long int* n, const long int* nrhs, NT2_WRAP_COMPLEX* ap, long int* ipiv, NT2_WRAP_COMPLEX* b, const long int* ldb, long int* info);
 void NT2_F77NAME(cpbsv)(const char* uplo, const long int* n, const long int* kd, const long int* nrhs, NT2_WRAP_COMPLEX* ab, const long int* ldab, NT2_WRAP_COMPLEX* b, const long int* ldb, long int* info);
 void NT2_F77NAME(cposv)(const char* uplo, const long int* n, const long int* nrhs, NT2_WRAP_COMPLEX* a, const long int* lda, NT2_WRAP_COMPLEX* b, const long int* ldb, long int* info);
 void NT2_F77NAME(cppsv)(const char* uplo, const long int* n, const long int* nrhs, NT2_WRAP_COMPLEX* ap, NT2_WRAP_COMPLEX* b, const long int* ldb, long int* info);
 void NT2_F77NAME(cptsv)(const long int* n, const long int* nrhs, float* d, NT2_WRAP_COMPLEX* e, NT2_WRAP_COMPLEX* b, const long int* ldb, long int* info);
 void NT2_F77NAME(cspsv)(const char* uplo, const long int* n, const long int* nrhs, NT2_WRAP_COMPLEX* ap, long int* ipiv, NT2_WRAP_COMPLEX* b, const long int* ldb, long int* info);
 void NT2_F77NAME(csysv)(const char* uplo, const long int* n, const long int* nrhs, NT2_WRAP_COMPLEX* a, const long int* lda, long int* ipiv, NT2_WRAP_COMPLEX* b, const long int* ldb, NT2_WRAP_COMPLEX* work, const long int* lwork, long int* info);

 void NT2_F77NAME(dgbsv)(const long int* n, const long int* kl, const long int* ku, const long int* nrhs, double* ab, const long int* ldab, long int* ipiv, double* b, const long int* ldb, long int* info);
 void NT2_F77NAME(dgesv)(const long int* n, const long int* nrhs, double* a, const long int* lda, long int* ipiv, double* b, const long int* ldb, long int* info);
 void NT2_F77NAME(dgtsv)(const long int* n, const long int* nrhs, double* dl, double* d, double* du, double* b, const long int* ldb, long int* info);
 void NT2_F77NAME(dpbsv)(const char* uplo, const long int* n, const long int* kd, const long int* nrhs, double* ab, const long int* ldab, double* b, const long int* ldb, long int* info);
 void NT2_F77NAME(dposv)(const char* uplo, const long int* n, const long int* nrhs, double* a, const long int* lda, double* b, const long int* ldb, long int* info);
 void NT2_F77NAME(dppsv)(const char* uplo, const long int* n, const long int* nrhs, double* ap, double* b, const long int* ldb, long int* info);
 void NT2_F77NAME(dptsv)(const long int* n, const long int* nrhs, double* d, double* e, double* b, const long int* ldb, long int* info);
 void NT2_F77NAME(dspsv)(const char* uplo, const long int* n, const long int* nrhs, double* ap, long int* ipiv, double* b, const long int* ldb, long int* info);
 void NT2_F77NAME(dsysv)(const char* uplo, const long int* n, const long int* nrhs, double* a, const long int* lda, long int* ipiv, double* b, const long int* ldb, double* work, const long int* lwork, long int* info);

 void NT2_F77NAME(sgbsv)(const long int* n, const long int* kl, const long int* ku, const long int* nrhs, float* ab, const long int* ldab, long int* ipiv, float* b, const long int* ldb, long int* info);
 void NT2_F77NAME(sgesv)(const long int* n, const long int* nrhs, float* a, const long int* lda, long int* ipiv, float* b, const long int* ldb, long int* info);
 void NT2_F77NAME(sgtsv)(const long int* n, const long int* nrhs, float* dl, float* d, float* du, float* b, const long int* ldb, long int* info);
 void NT2_F77NAME(spbsv)(const char* uplo, const long int* n, const long int* kd, const long int* nrhs, float* ab, const long int* ldab, float* b, const long int* ldb, long int* info);
 void NT2_F77NAME(sposv)(const char* uplo, const long int* n, const long int* nrhs, float* a, const long int* lda, float* b, const long int* ldb, long int* info);
 void NT2_F77NAME(sppsv)(const char* uplo, const long int* n, const long int* nrhs, float* ap, float* b, const long int* ldb, long int* info);
 void NT2_F77NAME(sptsv)(const long int* n, const long int* nrhs, float* d, float* e, float* b, const long int* ldb, long int* info);
 void NT2_F77NAME(sspsv)(const char* uplo, const long int* n, const long int* nrhs, float* ap, long int* ipiv, float* b, const long int* ldb, long int* info);
 void NT2_F77NAME(ssysv)(const char* uplo, const long int* n, const long int* nrhs, float* a, const long int* lda, long int* ipiv, float* b, const long int* ldb, float* work, const long int* lwork, long int* info);

 void NT2_F77NAME(zgbsv)(const long int* n, const long int* kl, const long int* ku, const long int* nrhs, NT2_WRAP_COMPLEX* ab, const long int* ldab, long int* ipiv, NT2_WRAP_COMPLEX* b, const long int* ldb, long int* info);
 void NT2_F77NAME(zgesv)(const long int* n, const long int* nrhs, NT2_WRAP_COMPLEX* a, const long int* lda, long int* ipiv, NT2_WRAP_COMPLEX* b, const long int* ldb, long int* info);
 void NT2_F77NAME(zgtsv)(const long int* n, const long int* nrhs, NT2_WRAP_COMPLEX* dl, NT2_WRAP_COMPLEX* d, NT2_WRAP_COMPLEX* du, NT2_WRAP_COMPLEX* b, const long int* ldb, long int* info);
 void NT2_F77NAME(zhesv)(const char* uplo, const long int* n, const long int* nrhs, NT2_WRAP_COMPLEX* a, const long int* lda, long int* ipiv, NT2_WRAP_COMPLEX* b, const long int* ldb, NT2_WRAP_COMPLEX* work, const long int* lwork, long int* info);
 void NT2_F77NAME(zhpsv)(const char* uplo, const long int* n, const long int* nrhs, NT2_WRAP_COMPLEX* ap, long int* ipiv, NT2_WRAP_COMPLEX* b, const long int* ldb, long int* info);
 void NT2_F77NAME(zpbsv)(const char* uplo, const long int* n, const long int* kd, const long int* nrhs, NT2_WRAP_COMPLEX* ab, const long int* ldab, NT2_WRAP_COMPLEX* b, const long int* ldb, long int* info);
 void NT2_F77NAME(zposv)(const char* uplo, const long int* n, const long int* nrhs, NT2_WRAP_COMPLEX* a, const long int* lda, NT2_WRAP_COMPLEX* b, const long int* ldb, long int* info);
 void NT2_F77NAME(zppsv)(const char* uplo, const long int* n, const long int* nrhs, NT2_WRAP_COMPLEX* ap, NT2_WRAP_COMPLEX* b, const long int* ldb, long int* info);
 void NT2_F77NAME(zptsv)(const long int* n, const long int* nrhs, double* d, NT2_WRAP_COMPLEX* e, NT2_WRAP_COMPLEX* b, const long int* ldb, long int* info);
 void NT2_F77NAME(zspsv)(const char* uplo, const long int* n, const long int* nrhs, NT2_WRAP_COMPLEX* ap, long int* ipiv, NT2_WRAP_COMPLEX* b, const long int* ldb, long int* info);
 void NT2_F77NAME(zsysv)(const char* uplo, const long int* n, const long int* nrhs, NT2_WRAP_COMPLEX* a, const long int* lda, long int* ipiv, NT2_WRAP_COMPLEX* b, const long int* ldb, NT2_WRAP_COMPLEX* work, const long int* lwork, long int* info);
  
}


#endif

// /////////////////////////////////////////////////////////////////////////////
// End of l_sv_C.hpp
// /////////////////////////////////////////////////////////////////////////////
