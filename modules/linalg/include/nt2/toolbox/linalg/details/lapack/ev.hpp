/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_DETAILS_LAPACK_EV_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_DETAILS_LAPACK_EV_HPP_INCLUDED
#include <nt2/toolbox/linalg/details/utility/f77_wrapper.hpp>
#include <nt2/toolbox/linalg/details/lapack/workspace.hpp>
#include <nt2/include/functions/max.hpp>
#include <nt2/table.hpp>
// syev,  heev

extern "C"
{
  void NT2_F77NAME(dsyev)(const char* jobz, const char* uplo, const la_int* n,
                          double* a, const la_int* lda, double* w, double* work,
                          const la_int* lwork, la_int* info);
  void NT2_F77NAME(ssyev)(const char* jobz, const char* uplo, const la_int* n,
                          float* a, const la_int* lda, float* w, float* work,
                          const la_int* lwork, la_int* info);
  void NT2_F77NAME(zheev)(const char* jobz, const char* uplo, const la_int* n,
                          la_complex* a, const la_int* lda, double* w, la_complex* work,
                          const la_int* lwork, double* rwork, la_int* info);
  void NT2_F77NAME(cheev)(const char* jobz, const char* uplo, const la_int* n,
                          la_complex* a, const la_int* lda, float* w, la_complex* work,
                          const la_int* lwork, float* rwork, la_int* info);
}

namespace nt2
{
  namespace details
  {
#define NT2_EV(NAME, T)                                                 \
    inline void ev(const char *jobz, const char*uplo, const la_int*n,   \
                   T* S, const la_int *lda,                             \
                   T *eigvals,                                          \
                   la_int *info)                                        \
    {                                                                   \
      table<T> w( nt2::of_size(1, 1) );                                 \
      la_int ldw = -1;                                                  \
      NT2_F77NAME(NAME)(jobz, uplo, n, s, lda, eigvals,                 \
                        w.raw(), &ldw, info);                           \
      ldw = (la_int)(w.raw[0]);                                         \
      w.resize(ofSize(ldw, 1));                                         \
      ldw =  nt2::leading_size(w);                                      \
      NT2_F77NAME(NAME)(jobz, uplo, n, s, lda, eigvals,                 \
                        w.raw(), &ldw, info);                           \
    }                                                                   \
     
    NT2_EV(dsyev, double)
    NT2_EV(ssyev, float)
      
#undef NT2_EV      
      
#define NT2_EV(NAME, T, TBASE)                                          \
    inline void ev(const char *jobz, const char*uplo, const la_int*n,   \
                   T* S, const la_int *lda,                             \
                   TBASE  *eigvals,                                     \
                   la_int *info)                                        \
    {                                                                   \
      table<T> w( nt2::of_size(1, 1) );                                 \
      la_int ldw = -1;                                                  \
      table<TBASE>  rw(nt2::of_size(nt2::max(1l, 3*(*n)-2), 1l));       \
      NT2_F77NAME(NAME)(jobz, uplo, N, S, lda, eigvals, w.begin()       \
                        , &ldw, rw.raw(), info);                        \
      ldw = (la_int)(w.raw[0]);                                         \
      w.resize(ofSize(ldw, 1l));                                        \
      ldw =  nt2::leading_size(w);                                      \
      NT2_F77NAME(NAME)(jobz, uplo, N, S, lda, eigvals, w.begin(),      \
                        &ldw, rw.raw(), info);                          \
    }                                                                   \
        
    NT2_EV(zheev, std::complex < double >, double)
    NT2_EV(cheev, std::complex < float >, float)  
      
#undef NT2_EV      

      }
}


#endif

