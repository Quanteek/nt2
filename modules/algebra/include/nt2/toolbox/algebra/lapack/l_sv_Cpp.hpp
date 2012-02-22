/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_ALGEBRA_LAPACK_L_SV_CPP_HPP_INCLUDED
#define NT2_TOOLBOX_ALGEBRA_LAPACK_L_SV_CPP_HPP_INCLUDED

#define NT2_WRAP_COMPLEX void
  
#define NT2_SV(T, PREFIX)                                               \
  void gbsv(const long int* n, const long int* kl, const long int* ku,  \
            const long int* nrhs,                                       \
            T* ab, const long int* ldab,                                \
            long int* ipiv, T* b, const long int* ldb,                  \
            long int* info);                                            \
  {                                                                     \
    BOOST_PP_CAT(PREFIX,NT2_F77NAME(gbsv,_))                            \
      (n,kl,ku,nrhs,ab,ldab,ipiv,b,ldb,info);                           \
  }                                                                     \
  /**/
  NT2_SV(double, d) 
  NT2_SV(float,  s) 
  NT2_SV(std::complex<double>, z) 
  NT2_SV(std::complex<float>, c) 
#undef NT2_SV

#define NT2_SV(T, PREFIX)                              \
  void gbsv(const long int* n, const long int* nrhs,   \
            T* a, const long int* lda,                 \
            long int* ipiv,                            \
            T* b, const long int* ldb,                 \
            long int* info);                           \
  {                                                    \
    BOOST_PP_CAT(PREFIX,NT2_F77NAME(gesv,_))           \
      (n,nrhs,ab,ldab,ipiv,b,ldb,info);                \
  }                                                    \
  /**/
  NT2_SV(double, d) 
  NT2_SV(float,  s) 
  NT2_SV(std::complex<double>, z) 
  NT2_SV(std::complex<float>, c) 
#undef NT2_SV
    
#define NT2_SV(T, PREFIX)                              \
  void gtsv(const long int* n, const long int* nrhs,   \
            T* dl, T* d, T* du,                        \
            T* b, const long int* ldb,                 \
            long int* info);                           \
  {                                                    \
    BOOST_PP_CAT(PREFIX,NT2_F77NAME(gtsv,_))           \
      (n,nrhs,dl,d,du,b,ldb,info);                     \
  }                                                    \
  /**/
  NT2_SV(double, d) 
  NT2_SV(float,  s) 
  NT2_SV(std::complex<double>, z) 
  NT2_SV(std::complex<float>, c) 
#undef NT2_SV
    
#define NT2_SV(T, PREFIX)                              \
  void hesv(const char* uplo,                          \
            const long int* n,                         \
            const long int* nrhs,                      \
            T* a, const long int* lda,                 \
            long int* ipiv,                            \
            T* b, const long int* ldb,                 \
            T* work, const long int* lwork,            \
            long int* info);                           \
  {                                                    \
    BOOST_PP_CAT(PREFIX,NT2_F77NAME(hesv,_))           \
      (uplo,n,nrhs,ab,ldab,ipiv,b,ldb,work,lwork,info); \
  }                                                     \
  /**/
  NT2_SV(double, d) 
  NT2_SV(float,  s) 
  NT2_SV(std::complex<double>, z) 
  NT2_SV(std::complex<float>, c) 
#undef NT2_SV
        
#define NT2_SV(T, PREFIX)                              \
  void hpsv(const char* uplo,                          \
            const long int* n,                         \
            const long int* nrhs,                      \
            T* a, const long int* lda,                 \
            long int* ipiv,                            \
            T* b, const long int* ldb,                 \
            long int* info);                           \
  {                                                    \
    BOOST_PP_CAT(PREFIX,NT2_F77NAME(hpsv,_))           \
      (uplo,n,nrhs,ab,ldab,ipiv,b,ldb,info);           \
  }                                                    \
  /**/
  NT2_SV(double, d) 
  NT2_SV(float,  s) 
  NT2_SV(std::complex<double>, z) 
  NT2_SV(std::complex<float>, c) 
#undef NT2_SV


#define NT2_SV(T, PREFIX)                              \
  void pbsv(const char* uplo,                          \
            const long int* n,                         \
            const long int* kd,                        \
            const long int* nrhs,                      \
            T* a, const long int* lda,                 \
            long int* ipiv,                            \
            T* b, const long int* ldb,                 \
            long int* info);                           \
  {                                                    \
    BOOST_PP_CAT(PREFIX,NT2_F77NAME(pbsv,_))           \
      (uplo,n,kdnrhs,ab,ldab,ipiv,b,ldb,info);         \
  }                                                    \
  /**/
  NT2_SV(double, d) 
  NT2_SV(float,  s) 
  NT2_SV(std::complex<double>, z) 
  NT2_SV(std::complex<float>, c) 
#undef NT2_SV

#define NT2_SV(T, PREFIX)                              \
  void posv(const char* uplo,                          \
            const long int* n,                         \
            const long int* nrhs,                      \
            T* a, const long int* lda,                 \
            long int* ipiv,                            \
            T* b, const long int* ldb,                 \
            long int* info);                           \
  {                                                    \
    BOOST_PP_CAT(PREFIX,NT2_F77NAME(posv,_))           \
      (uplo,n,nrhs,ab,ldab,ipiv,b,ldb,info);            \
  }                                                     \
  /**/
  NT2_SV(double, d) 
  NT2_SV(float,  s) 
  NT2_SV(std::complex<double>, z) 
  NT2_SV(std::complex<float>, c) 
#undef NT2_SV
        
#define NT2_SV(T, PREFIX)                              \
  void ppsv(const char* uplo,                          \
            const long int* n,                         \
            const long int* nrhs,                      \
            T* ap,                                     \
            long int* ipiv,                            \
            T* b, const long int* ldb,                 \
            long int* info);                           \
  {                                                    \
    BOOST_PP_CAT(PREFIX,NT2_F77NAME(ppsv,_))           \
      (uplo,n,nrhs,ab,ldab,ipiv,b,ldb,info);            \
  }                                                     \
  /**/
  NT2_SV(double, d) 
  NT2_SV(float,  s) 
  NT2_SV(std::complex<double>, z) 
  NT2_SV(std::complex<float>, c) 
#undef NT2_SV
    
#define NT2_SV(T, BASET, PREFIX)                       \
  void ptsv(const char* uplo,                          \
            const long int* n,                         \
            BASET* d,                                  \
            T* e,                                      \
            long int* ipiv,                            \
            T* b, const long int* ldb,                 \
            long int* info);                           \
  {                                                    \
    BOOST_PP_CAT(PREFIX,NT2_F77NAME(ptsv,_))           \
      (uplo,n,nrhs,ab,ldab,ipiv,b,ldb,info);            \
  }                                                     \
  /**/
  NT2_SV(double, double, d) 
  NT2_SV(float, float, s) 
  NT2_SV(std::complex<double>, double, z) 
  NT2_SV(std::complex<float>,float, c) 
#undef NT2_SV
        
#define NT2_SV(T, PREFIX)                              \
  void spsv(const char* uplo,                          \
            const long int* n,                         \
            const long int* nrhs,                      \
            T* ap,                                     \
            long int* ipiv,                            \
            T* b, const long int* ldb,                 \
            long int* info);                           \
  {                                                    \
    BOOST_PP_CAT(PREFIX,NT2_F77NAME(spsv,_))           \
      (uplo,n,nrhs,ab,ldab,ipiv,b,ldb,info);           \
  }                                                    \
  /**/
  NT2_SV(double, d) 
  NT2_SV(float,  s) 
  NT2_SV(std::complex<double>, z) 
  NT2_SV(std::complex<float>, c) 
#undef NT2_SV

#define NT2_SV(T, BASET, PREFIX)                       \
  void sysv(const char* uplo,                          \
            const long int* n,                         \
            const long int* nrhs,                      \
            T* a, const long int* lda,                 \
            long int* ipiv,                            \
            T* b, const long int* ldb,                 \
            T* work, const long int* lwork,            \
            long int* info);                           \
  {                                                    \
    BOOST_PP_CAT(PREFIX,NT2_F77NAME(sysv,_))           \
      (uplo,n,nrhs,ab,ldab,ipiv,b,ldb,work,lwork,info); \
  }                                                     \
  /**/
  NT2_SV(double, double, d) 
  NT2_SV(float, float, s) 
  NT2_SV(std::complex<double>, double, z) 
  NT2_SV(std::complex<float>,float, c) 
#undef NT2_SV

}


#endif

// /////////////////////////////////////////////////////////////////////////////
// End of l_sv_C.hpp
// /////////////////////////////////////////////////////////////////////////////
