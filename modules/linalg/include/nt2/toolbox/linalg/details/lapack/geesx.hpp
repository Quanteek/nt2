/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_DETAILS_LAPACK_GEESX_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_DETAILS_LAPACK_GEESX_HPP_INCLUDED
// tresx

namespace nt2
{
  namespace details
  {

    inline long int selectall(void* ) { return true; }
    inline long int selectall2( const void*, const  void* ) { return true; }

    
    extern "C"
    {
#define NT2_COMPLEX void
      void F77NAME(cgeesx)(const char* jobvs, const char* sort, const selectall* select, const char* sense, const long int* n, COMPLEX* a, const long int* lda, long int* sdim, COMPLEX* w, const COMPLEX* vs, const long int* ldvs, float* rconde, float* rcondv, COMPLEX* work, const long int* lwork, float* rwork, long int* bwork, long int* info);
      void F77NAME(sgeesx)(const char* jobvs, const char* sort, const selectall2* selectall , const char* sense, const long int* n, float* a, const long int* lda, long int* sdim, float* wr, float* wi, const float* vs, const long int* ldvs, float* rconde, float* rcondv, float* work, const long int* lwork, long int* iwork, const long int* liwork, long int* bwork, long int* info);
      void F77NAME(zgeesx)(const char* jobvs, const char* sort, const lselectall* select , const char* sense, const long int* n, COMPLEX* a, const long int* lda, long int* sdim, COMPLEX* w, const COMPLEX* vs, const long int* ldvs, double* rconde, double* rcondv, COMPLEX* work, const long int* lwork, double* rwork, long int* bwork, long int* info);
      void F77NAME(dgeesx)(const char* jobvs, const char* sort, const selectall2*select , const char* sense, const long int* n, double* a, const long int* lda, long int* sdim, double* wr, double* wi, const double* vs, const long int* ldvs, double* rconde, double* rcondv, double* work, const long int* lwork, long int* iwork, const long int* liwork, long int* bwork, long int* info);
#undef NT2_COMPLEX
    }

#define NT2_GEESX(NAME, T, TBASE)               \
    inline void geesx(const char* jobvs,        \
                      const char* sort,         \
                      const long int* select,   \
                      const char* sense,        \
                      const long int* n,        \
                      T* a,                     \
                      const long int* lda,      \
                      long int* sdim,           \
                      T* ws,                    \
                      const T* vs,              \
                      const long int* ldvs,     \
                      TBASE* rconde,            \
                      TBASE* rcondv,            \
                      long int* info,           \
                      workspace<T> & w)         \
    {                                           \
      w.resizerw(*n);                                            \
      w.resizebw((lower(*sort) == 'n') ? 0 : *n);                \
      NT2_F77NAME( NAME )(jobvs, sort, select, sense, n,         \
                          a, lda, sdim, ws, vs, ldvs,            \
                          rconde, rcondv, w.getw(), w.query(),   \
                          w.getrw(), w.getbw(), info);           \
      w.resizew(w.neededsize());                                 \
      NT2_F77NAME( NAME )(jobvs, sort, select, sense, n,         \
                          a, lda, sdim, ws, vs, ldvs,            \
                          rconde, rcondv, w.getw(),              \
                          &w.neededsize(), w.getrw(),            \
                          w.getbw(), info);                      \
    }                                                            \
    inline void geesx(const char* jobvs,                         \
                      const char* sort,                          \
                      const long int* select,                    \
                      const char* sense,                         \
                      const long int* n,                         \
                      T* a,                                      \
                      const long int* lda,                       \
                      long int* sdim,                            \
                      T* ws,                                     \
                      const T* vs,                               \
                      const long int* ldvs,                      \
                      TBASE* rconde,                             \
                      TBASE* rcondv,                             \
                      long int* info)                            \
    {                                                            \
      workspace<T> w;                                            \
      geesx(jobvs, sort, select, sense, n,                       \
            a, lda, sdim, ws, vs, ldvs, rconde,                  \
            rcondv, info, w);                                    \
    }                                                            \

    NT2_GEESX(sgeesx, float,  float)
    NT2_GEESX(dgeesx, double, double)
    NT2_GEESX(cgeesx, std::complex<float>,  float)
    NT2_GEESX(zgeesx, std::complex<double>, double)

#undef NT2_GEESX


}

#endif
