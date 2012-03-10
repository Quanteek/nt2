/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_DETAILS_LAPACK_GEEVX_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_DETAILS_LAPACK_GEEVX_HPP_INCLUDED
#include <nt2/toolbox/linalg/details/utility/f77_wrapper.hpp>
#include <nt2/toolbox/linalg/details/lapack/workspace.hpp>
// gecon,  geevx


namespace nt2
{
  namespace details
  {
    //////////////////////////////////////////////////////////////////////
    // geevx calls
    //////////////////////////////////////////////////////////////////////

    extern "C"
    {
#define NT2_COMPLEX void
      void NT2_F77NAME(cgeevx)(const char* balanc, const char* jobvl, const char* jobvr,
                           const char* sense, const long int* n, const
                           NT2_COMPLEX* a, const long int* lda, NT2_COMPLEX* w,
                           const NT2_COMPLEX* vl, const long int* ldvl,
                           const NT2_COMPLEX* vr, const long int* ldvr, long int* ilo, long int* ihi,
                           float* scale, float* abnrm, float* rconde, float* rcondv,
                           NT2_COMPLEX* work, const long int* lwork, float* rwork, long int* info);
      void NT2_F77NAME(dgeevx)(const char* balanc, const char* jobvl, const char* jobvr,
                           const char* sense, const long int* n,
                           const double* a, const long int* lda, double* wr, double* wi,
                           const double* vl, const long int* ldvl,
                           const double* vr, const long int* ldvr, long int* ilo, long int* ihi,
                           double* scale, double* abnrm, double* rconde, double* rcondv,
                           double* work, const long int* lwork, long int* iwork, long int* info);
      void NT2_F77NAME(sgeevx)(const char* balanc, const char* jobvl, const char* jobvr,
                           const char* sense, const long int* n,
                           const float* a, const long int* lda, float* wr, float* wi,
                           const float* vl, const long int* ldvl,
                           const float* vr, const long int* ldvr,
                           long int* ilo, long int* ihi, float* scale, float* abnrm,
                           float* rconde, float* rcondv,
                           float* work, const long int* lwork, long int* iwork, long int* info);
      void NT2_F77NAME(zgeevx)(const char* balanc, const char* jobvl, const char* jobvr,
                           const char* sense, const long int* n,
                           const NT2_COMPLEX* a, const long int* lda, NT2_COMPLEX* w,
                           const NT2_COMPLEX* vl, const long int* ldvl,
                           const NT2_COMPLEX* vr, const long int* ldvr,
                           long int* ilo, long int* ihi, double* scale, double* abnrm,
                           double* rconde, double* rcondv,
                           NT2_COMPLEX* work, const long int* lwork, double* rwork, long int* info);
      #undef NT2_COMPLEX
    }

#define NT2_GEEVX(NAME, T)                                      \
    inline void geevx(const char* balanc,                       \
                      const char* jobvl,                        \
                      const char* jobvr,                        \
                      const char* sense,                        \
                      const long int* n,                        \
                      const T* a,                               \
                      const long int* lda,                      \
                      T* wr,                                    \
                      T* wi,                                    \
                      const T* vl,                              \
                      const long int* ldvl,                     \
                      const T* vr,                              \
                      const long int* ldvr,                     \
                      long int* ilo,                            \
                      long int* ihi,                            \
                      T* scale,                                 \
                      T* abnrm,                                 \
                      T* rconde,                                \
                      T* rcondv,                                \
                      long int* info,                           \
                      nt2::details::workspace<T> & w)           \
    {                                                           \
      w.resizeiw((nt2::details::lower(*sense) == 'n' ||         \
                  nt2::details::lower(*sense) ==  'e')          \
                 ? 0                                            \
                 : 2**n-2                                       \
                 );                                             \
      NT2_F77NAME( NAME )(balanc, jobvl, jobvr, sense,          \
                          n, a, lda, wr, wi, vl, ldvl,          \
                          vr, ldvr, ilo, ihi, scale,            \
                          abnrm, rconde, rcondv, w.getw(),      \
                          w.query(), w.getiw(), info);          \
      w.resizew(w.neededsize());                                \
      NT2_F77NAME( NAME )(balanc, jobvl, jobvr, sense,          \
                          n, a, lda, wr, wi, vl, ldvl,          \
                          vr, ldvr, ilo, ihi, scale,            \
                          abnrm, rconde, rcondv, w.getw(),      \
                          &w.neededsize(), w.getiw(),           \
                          info);                                \
    }                                                           \
    inline void geevx(const char* balanc,                       \
                      const char* jobvl,                        \
                      const char* jobvr,                        \
                      const char* sense,                        \
                      const long int* n,                        \
                      const T* a,                               \
                      const long int* lda,                      \
                      T* wr,                                    \
                      T* wi,                                    \
                      const T* vl,                              \
                      const long int* ldvl,                     \
                      const T* vr,                              \
                      const long int* ldvr,                     \
                      long int* ilo,                            \
                      long int* ihi,                            \
                      T* scale,                                 \
                      T* abnrm,                                 \
                      T* rconde,                                \
                      T* rcondv,                                \
                      long int* info)                           \
    {                                                           \
      nt2::details::workspace<T> w;                             \
      geevx(balanc, jobvl, jobvr, sense,                        \
            n, a, lda, wr, wi, vl, ldvl,                        \
            vr, ldvr, ilo, ihi, scale,                          \
            abnrm, rconde, rcondv,                              \
            info, w);                                           \
    }                                                           \

    NT2_GEEVX(sgeevx, float)
    NT2_GEEVX(dgeevx, double)

#undef NT2_GEEVX
      
#define NT2_GEEVX(NAME, T, TBASE)                               \
      inline void geevx(const char* balanc,                     \
                        const char* jobvl,                      \
                        const char* jobvr,                      \
                        const char* sense,                      \
                        const long int* n,                      \
                        const T* a,                             \
                        const long int* lda,                    \
                        T* ws,                                  \
                        const T* vl,                            \
                        const long int* ldvl,                   \
                        const T* vr,                            \
                        const long int* ldvr,                   \
                        long int* ilo,                          \
                        long int* ihi,                          \
                        TBASE* scale,                           \
                        TBASE* abnrm,                           \
                        TBASE* rconde,                          \
                        TBASE* rcondv,                          \
                        long int* info,                         \
                        nt2::details::workspace<T> & w)         \
      {                                                         \
        w.resizerw(2**n);                                       \
        F77NAME( NAME )(balanc, jobvl, jobvr, sense,            \
                        n, a, lda, ws, vl, ldvl, vr, ldvr,      \
                        ilo, ihi, scale, abnrm, rconde,         \
                        rcondv, w.getw(), w.query(), w.getrw(), \
                        info);                                  \
        w.resizew(w.neededsize());                              \
        F77NAME( NAME )(balanc, jobvl, jobvr, sense,            \
                        n, a, lda, ws, vl, ldvl, vr, ldvr,      \
                        ilo, ihi, scale, abnrm, rconde,         \
                        rcondv, w.getw(),                       \
                        &w.neededsize(), w.getrw(),             \
                        info);                                  \
      }                                                         \
          inline void geevx(                                    \
                            const char* balanc,                 \
                            const char* jobvl,                  \
                            const char* jobvr,                  \
                            const char* sense,                  \
                            const long int* n,                  \
                            const T* a,                         \
                            const long int* lda,                \
                            T* ws,                              \
                            const T* vl,                        \
                            const long int* ldvl,               \
                            const T* vr,                        \
                            const long int* ldvr,               \
                            long int* ilo,                      \
                            long int* ihi,                      \
                            TBASE* scale,                       \
                            TBASE* abnrm,                       \
                            TBASE* rconde,                      \
                            TBASE* rcondv,                      \
                            long int* info)                     \
          {                                                     \
            nt2::details::workspace<T> w;                       \
            geevx(balanc, jobvl, jobvr, sense,                  \
                  n, a, lda, ws, vl, ldvl, vr, ldvr,            \
                  ilo, ihi, scale, abnrm, rconde,               \
                  rcondv, info, w);                             \
          }                                                     \
          
    NT2_GEEVX(cgeevx, std::complex<float>, float)
    NT2_GEEVX(zgeevx, std::complex<double>, double)

#undef NT2_GEEVX
      }
}


#endif

// /////////////////////////////////////////////////////////////////////////////
// End of con.hpp
// /////////////////////////////////////////////////////////////////////////////
