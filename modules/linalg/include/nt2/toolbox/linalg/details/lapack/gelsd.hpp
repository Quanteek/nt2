/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_DETAILS_LAPACK_GELSD_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_DETAILS_LAPACK_GELSD_HPP_INCLUDED
#include <nt2/toolbox/linalg/details/utility/f77_wrapper.hpp>
#include <nt2/toolbox/linalg/details/lapack/utility.hpp>

namespace nt2
{
  namespace details
  {
    extern "C"
    {
#define NT2_COMPLEX void
  void NT2_F77NAME(cgelsd)(const long int* m, const long int* n, const long int* nrhs,
                           NT2_COMPLEX* a, const long int* lda, NT2_COMPLEX* b, const long int* ldb,
                           const float* s, const float* rcond, long int* rank, NT2_COMPLEX* work,
                           const long int* lwork, float* rwork, long int* iwork, long int* info);
  void NT2_F77NAME(dgelsd)(const long int* m, const long int* n, const long int* nrhs,
                           double* a, const long int* lda, double* b, const long int* ldb,
                           const double* s, const double* rcond, long int* rank, double* work,
                           const long int* lwork, long int* iwork, long int* info);
  void NT2_F77NAME(sgelsd)(const long int* m, const long int* n, const long int* nrhs,
                           float* a, const long int* lda, float* b, const long int* ldb,
                           const float* s, const float* rcond, long int* rank, float* work,
                           const long int* lwork, long int* iwork, long int* info);
  void NT2_F77NAME(zgelsd)(const long int* m, const long int* n, const long int* nrhs,
                           NT2_COMPLEX* a, const long int* lda, NT2_COMPLEX* b, const long int* ldb,
                           const double* s, const double* rcond, long int* rank, NT2_COMPLEX* work,
                           const long int* lwork, double* rwork, long int* iwork, long int* info);
#undef NT2_COMPLEX
    }

    struct gelsdUtils
    {
      static inline long int nlvl(const std::string & name, long int minnm, long int nrsh)
      {
        std::string blank =  " ";
        minnm = std::max(minnm, 1l);
        long int smlsiz = nt2::details::EnvBlockSize(9, name.c_str(), blank,0, 0, 0, 0);
        return (long int)nt2::max( int(nt2::log(double(minnm))/ double(smlsiz+1))/nt2::log(2.0)+1,0.0);
      }      
    
      static inline size_t liwork(const std::string & name, long int n,  long int m, long int nrsh)
      {
        long int minnm = nt2::min(n, m); 
        return (3*nlvl(name, minnm, nrsh)+11)*minnm; 
      }      
      
      static inline long int lrwork(const std::string & name, long int n,  long int m, long int nrsh)
      {
        std::string blank =  " ";
        long int maxnm = nt2::max(n, m);
        long int minnm = nt2::min(n, m); 
        long int nl = nlvl(name, minnm, nrsh); 
        long int smlsiz = nt2::details::EnvBlockSize(9, name.c_str(), blank,0, 0, 0, 0);
        return 10*maxnm + 2*maxnm*smlsiz + 8*m*nl + 3*smlsiz*nrsh +  (smlsiz+1)*(smlsiz+1); 
      }      
    };    
    
    
#define NT2_GELSD(NAME, T)                      \
    inline void gelsd(const long int* m,        \
                      const long int* n,        \
                      const long int* nrhs,     \
                      T* a,                     \
                      const long int* lda,      \
                      T* b,                     \
                      const long int* ldb,      \
                      const T* s,               \
                      const T* rcond,           \
                      long int* rank,           \
                      long int* info,           \
              nt2::details::workspace<T> & w)   \
    {                                           \
      w.resizeiw(gelsdUtils::liwork("AME", *n, *m, *nrhs));     \
      NT2_F77NAME( NAME )(m, n, nrhs,           \
                  a, lda, b, ldb,               \
                  s, rcond, rank, w.getw(),     \
                  w.query(), w.getiw(),         \
                  info);                        \
      w.resizew(w.neededsize());                \
      NT2_F77NAME( NAME )(m, n, nrhs,           \
                  a, lda, b, ldb,               \
                  s, rcond, rank, w.getw(),     \
                  &w.neededsize(), w.getiw(),   \
                  info);                        \
    }                                           \
    inline void gelsd(const long int* m,        \
                      const long int* n,        \
                      const long int* nrhs,     \
                      T* a,                     \
                      const long int* lda,      \
                      T* b,                     \
                      const long int* ldb,      \
                      const T* s,               \
                      const T* rcond,           \
                      long int* rank,           \
                      long int* info)           \
    {                                           \
      nt2::details::workspace<T> w;             \
      gelsd(m, n, nrhs,                         \
            a, lda, b, ldb,                     \
            s, rcond, rank, info, w);           \
    }                                           \

    NT2_GELSD(sgelsd, float)
    NT2_GELSD(dgelsd, double)

#undef NT2_GELSD
      
#define NT2_GELSD(NAME, T, TBASE)               \
      inline void gelsd(const long int* m,      \
                        const long int* n,      \
                        const long int* nrhs,   \
                        T* a,                   \
                        const long int* lda,    \
                        T* b,                   \
                        const long int* ldb,    \
                        const TBASE* s,         \
                        const TBASE* rcond,     \
                        long int* rank,         \
                        long int* info,         \
                        nt2::details::workspace<T> & w) \
      {                                                               \
        w.resizeiw(gelsdUtils::liwork("NAME", *n, *m, *nrhs));        \
        w.resizerw(gelsdUtils::lrwork("NAME", *n, *m, *nrhs));        \
        NT2_F77NAME( NAME )(m, n, nrhs,                               \
                            a, lda, b, ldb,                           \
                            s, rcond, rank,                           \
                            w.getw(), w.query(), w.getrw(),           \
                            w.getiw(), info);                         \
        w.resizew(w.neededsize());                                    \
        NT2_F77NAME( NAME )(m, n, nrhs,                               \
                            a, lda, b, ldb,                           \
                            s, rcond, rank,                           \
                            w.getw(), &w.neededsize(),                \
                            w.getrw(), w.getiw(), info);              \
      }                                                               \
     inline void gelsd(const long int* m,                             \
                       const long int* n,                             \
                       const long int* nrhs,                          \
                       T* a,                                          \
                       const long int* lda,                           \
                       T* b,                                          \
                       const long int* ldb,                           \
                       const TBASE* s,                                \
                       const TBASE* rcond,                            \
                       long int* rank,                                \
                       long int* info)                                \
     {                                                                \
       nt2::details::workspace<T> w;                                  \
       gelsd(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, info, w);    \
     }                                                                \

    NT2_GELSD(cgelsd, std::complex<float>,  float)
    NT2_GELSD(zgelsd, std::complex<double>, double)

#undef NT2_GELSD


  }
}

#endif

