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

extern "C"
{
  void NT2_F77NAME(cgelsd)(const la_int* m, const la_int* n, const la_int* nrhs,
                           la_complex* a, const la_int* lda, la_complex* b, const la_int* ldb,
                           const float* s, const float* rcond, la_int* rank, la_complex* work,
                           const la_int* lwork, float* rwork, la_int* iwork, la_int* info);
  void NT2_F77NAME(dgelsd)(const la_int* m, const la_int* n, const la_int* nrhs,
                           double* a, const la_int* lda, double* b, const la_int* ldb,
                           const double* s, const double* rcond, la_int* rank, double* work,
                           const la_int* lwork, la_int* iwork, la_int* info);
  void NT2_F77NAME(sgelsd)(const la_int* m, const la_int* n, const la_int* nrhs,
                           float* a, const la_int* lda, float* b, const la_int* ldb,
                           const float* s, const float* rcond, la_int* rank, float* work,
                           const la_int* lwork, la_int* iwork, la_int* info);
  void NT2_F77NAME(zgelsd)(const la_int* m, const la_int* n, const la_int* nrhs,
                           la_complex* a, const la_int* lda, la_complex* b, const la_int* ldb,
                           const double* s, const double* rcond, la_int* rank, la_complex* work,
                           const la_int* lwork, double* rwork, la_int* iwork, la_int* info);
}
namespace nt2
{
  namespace details
  {


    struct gelsdUtils
    {
      static inline la_int nlvl(const std::string & name, la_int minnm, la_int nrsh)
      {
        std::string blank =  " ";
        minnm = std::max(minnm, 1l);
        la_int smlsiz = nt2::details::EnvBlockSize(9, name.c_str(), blank,0, 0, 0, 0);
        return (la_int)nt2::max( int(nt2::log(double(minnm))/ double(smlsiz+1))/nt2::log(2.0)+1,0.0);
      }      
    
      static inline size_t liwork(const std::string & name, la_int n,  la_int m, la_int nrsh)
      {
        la_int minnm = nt2::min(n, m); 
        return (3*nlvl(name, minnm, nrsh)+11)*minnm; 
      }      
      
      static inline la_int lrwork(const std::string & name, la_int n,  la_int m, la_int nrsh)
      {
        std::string blank =  " ";
        la_int maxnm = nt2::max(n, m);
        la_int minnm = nt2::min(n, m); 
        la_int nl = nlvl(name, minnm, nrsh); 
        la_int smlsiz = nt2::details::EnvBlockSize(9, name.c_str(), blank,0, 0, 0, 0);
        return 10*maxnm + 2*maxnm*smlsiz + 8*m*nl + 3*smlsiz*nrsh +  (smlsiz+1)*(smlsiz+1); 
      }      
    };    
    
    
#define NT2_GELSD(NAME, T)                      \
    inline void gelsd(const la_int* m,          \
                      const la_int* n,          \
                      const la_int* nrhs,       \
                      T* a,                     \
                      const la_int* lda,        \
                      T* b,                     \
                      const la_int* ldb,        \
                      const T* s,               \
                      const T* rcond,           \
                      la_int* rank,             \
                      la_int* info,             \
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
    inline void gelsd(const la_int* m,          \
                      const la_int* n,          \
                      const la_int* nrhs,       \
                      T* a,                     \
                      const la_int* lda,        \
                      T* b,                     \
                      const la_int* ldb,        \
                      const T* s,               \
                      const T* rcond,           \
                      la_int* rank,             \
                      la_int* info)             \
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
      inline void gelsd(const la_int* m,        \
                        const la_int* n,        \
                        const la_int* nrhs,     \
                        T* a,                   \
                        const la_int* lda,      \
                        T* b,                   \
                        const la_int* ldb,      \
                        const TBASE* s,         \
                        const TBASE* rcond,     \
                        la_int* rank,           \
                        la_int* info,           \
                nt2::details::workspace<T> & w) \
      {                                         \
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
      inline void gelsd(const la_int* m,                              \
                        const la_int* n,                              \
                        const la_int* nrhs,                           \
                        T* a,                                         \
                        const la_int* lda,                            \
                        T* b,                                         \
                        const la_int* ldb,                            \
                        const TBASE* s,                               \
                        const TBASE* rcond,                           \
                        la_int* rank,                                 \
                        la_int* info)                                 \
      {                                                               \
        nt2::details::workspace<T> w;                                 \
        gelsd(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, info, w);   \
      }                                                               \
          
    NT2_GELSD(cgelsd, std::complex<float>,  float)
    NT2_GELSD(zgelsd, std::complex<double>, double)

#undef NT2_GELSD


  }
}

#endif

