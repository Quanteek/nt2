/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_DETAILS_LAPACK_SVD_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_DETAILS_LAPACK_SVD_HPP_INCLUDED
#include <nt2/toolbox/linalg/details/utility/f77_wrapper.hpp>
#include <nt2/toolbox/linalg/details/lapack/workspace.hpp>

namespace nt2
{
  namespace details
  {
    /**  purpose
     **  =======
     **
     **  xgesvd computes the singular value decomposition (svd) of a DATA TYPE
     **  m-by-n matrix a, optionally computing the left and/or right singular
     **  vectors. the svd is written
     **
     **       a = u * sigma * conjugate-transpose(v)
     **
     **  where sigma is an m-by-n matrix which is zero except for its
     **  min(m,n) diagonal elements, u is an m-by-m unitary matrix, and
     **  v is an n-by-n unitary matrix.  the diagonal elements of sigma
     **  are the singular values of a; they are BASE DATA TYPE and non-negative, and
     **  are returned in descending order.  the first min(m,n) columns of
     **  u and v are the left and right singular vectors of a.
     **
     **  note that the routine returns v**h, not v.
     **
     **  arguments
     **  =========
     **
     **  jobu    (input) char
     **          specifies options for computing all or part of the matrix u:
     **          = 'a':  all m columns of u are returned in array u:
     **          = 's':  the first min(m,n) columns of u (the left singular
     **                  vectors) are returned in the array u;
     **          = 'o':  the first min(m,n) columns of u (the left singular
     **                  vectors) are overwritten on the array a;
     **          = 'n':  no columns of u (no left singular vectors) are
     **                  computed.
     **
     **  jobvt   (input) char
     **          specifies options for computing all or part of the matrix
     **          v**h:
     **          = 'a':  all n rows of v**h are returned in the array vt;
     **          = 's':  the first min(m,n) rows of v**h (the right singular
     **                  vectors) are returned in the array vt;
     **          = 'o':  the first min(m,n) rows of v**h (the right singular
     **                  vectors) are overwritten on the array a;
     **          = 'n':  no rows of v**h (no right singular vectors) are
     **                  computed.
     **
     **          jobvt and jobu cannot both be 'o'.
     **
     **  m       (input) long int
     **          the number of rows of the input matrix a.  m >= 0.
     **
     **  n       (input) long int
     **          the number of columns of the input matrix a.  n >= 0.
     **
     **  a       (input/output) DATA TYPE array, dimension (lda,n)
     **          on entry, the m-by-n matrix a.
     **          on exit,
     **          if jobu = 'o',  a is overwritten with the first min(m,n)
     **                          columns of u (the left singular vectors,
     **                          stored columnwise);
     **          if jobvt = 'o', a is overwritten with the first min(m,n)
     **                          rows of v**h (the right singular vectors,
     **                          stored rowwise);
     **          if jobu .ne. 'o' and jobvt .ne. 'o', the contents of a
     **                          are destroyed.
     **
     **  lda     (input) long int
     **          the leading dimension of the array a.  lda >= max(1,m).
     **
     **  s       (output) BASE DATA TYPE array, dimension (min(m,n))
     **          the singular values of a, sorted so that s(i) >= s(i+1).
     **
     **  u       (output) DATA TYPE array, dimension (ldu,ucol)
     **          (ldu,m) if jobu = 'a' or (ldu,min(m,n)) if jobu = 's'.
     **          if jobu = 'a', u contains the m-by-m unitary matrix u;
     **          if jobu = 's', u contains the first min(m,n) columns of u
     **          (the left singular vectors, stored columnwise);
     **          if jobu = 'n' or 'o', u is not referenced.
     **
     **  ldu     (input) long int
     **          the leading dimension of the array u.  ldu >= 1; if
     **          jobu = 's' or 'a', ldu >= m.
     **
     **  vt      (output) DATA TYPE array, dimension (ldvt,n)
     **          if jobvt = 'a', vt contains the n-by-n unitary matrix
     **          v**h;
     **          if jobvt = 's', vt contains the first min(m,n) rows of
     **          v**h (the right singular vectors, stored rowwise);
     **          if jobvt = 'n' or 'o', vt is not referenced.
     **
     **  ldvt    (input) long int
     **          the leading dimension of the array vt.  ldvt >= 1; if
     **          jobvt = 'a', ldvt >= n; if jobvt = 's', ldvt >= min(m,n).
     **
     **
     **
     **
     **  info    (output) long int
     **          = 0:  successful exit.
     **          < 0:  if info = -i, the i-th argument had an illegal value.
     **          > 0:  if cbdsqr did not converge, info specifies how many
     **                superdiagonals of an intermediate bidiagonal form b
     **                did not converge to zero. see the description of RWORK
     **                above for details.
     **
     **/
    extern "C"
    {
      #define NT2_COMPLEX void
      void NT2_F77NAME(cgesvd)(const char* jobu, const char* jobvt, const long int* m, const long int* n,
                           NT2_COMPLEX* a, const long int* lda, float* s, const NT2_COMPLEX* u, const long int* ldu,
                           const NT2_COMPLEX* vt, const long int* ldvt,
                           NT2_COMPLEX* work, const long int* lwork, float* rwork, long int* info);
      void NT2_F77NAME(dgesvd)(const char* jobu, const char* jobvt, const long int* m, const long int* n,
                           double* a, const long int* lda, double* s, const double* u, const long int* ldu,
                           const double* vt, const long int* ldvt,
                           double* work, const long int* lwork, long int* info);
      void NT2_F77NAME(sgesvd)(const char* jobu, const char* jobvt, const long int* m, const long int* n,
                           float* a, const long int* lda, float* s, const float* u, const long int* ldu,
                           const float* vt, const long int* ldvt,
                           float* work, const long int* lwork, long int* info);
      void NT2_F77NAME(zgesvd)(const char* jobu, const char* jobvt, const long int* m, const long int* n,
                           NT2_COMPLEX* a, const long int* lda, double* s, const NT2_COMPLEX* u, const long int* ldu,
                           const NT2_COMPLEX* vt, const long int* ldvt,
                           NT2_COMPLEX* work, const long int* lwork, double* rwork, long int* info);
      #undef NT2_COMPLEX
    }
    
#define NT2_GESVD(NAME, T, TBASE)                                       \
    inline void gesvd(const char* jobu,                                 \
                      const char* jobvt,                                \
                      const long int* m,                                \
                      const long int* n,                                \
                      T* a,                                             \
                      const long int* lda,                              \
                      TBASE* s,                                         \
                      const T* u,                                       \
                      const long int* ldu,                              \
                      const T* vt,                                      \
                      const long int* ldvt,                             \
                      long int* info,                                   \
                      nt2::details::workspace<T> & w)                   \
    {                                                                   \
      w.resizerw(5*nt2::min(*m,*n));                                    \
      NT2_F77NAME( NAME )(jobu, jobvt, m, n, a, lda, s, u, ldu,         \
                      vt, ldvt, w.getw(), w.query(), w.getrw(), info);  \
      w.resizew(w.neededsize());                                        \
      NT2_F77NAME( NAME )(jobu, jobvt, m, n, a, lda, s, u, ldu,         \
                      vt,ldvt,w.getw(),&w.neededsize(),w.getrw(),info); \
    }                                                                   \
    inline void gesvd(const char* jobu,                                 \
                      const char* jobvt,                                \
                      const long int* m,                                \
                      const long int* n,                                \
                      T* a,                                             \
                      const long int* lda,                              \
                      TBASE* s,                                         \
                      const T* u,                                       \
                      const long int* ldu,                              \
                      const T* vt,                                      \
                      const long int* ldvt,                             \
                      long int* info)                                   \
    {                                                                   \
      nt2::details::workspace<T> w;                                     \
      gesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, info, w);   \
    }                                                                   \

    NT2_GESVD(cgesvd, std::complex<float>,  float)
    NT2_GESVD(zgesvd, std::complex<double>, double)
#undef NT2_GESVD

#define NT2_GESVD(NAME, T)                                              \
    inline void gesvd(const char* jobu,                                 \
                      const char* jobvt,                                \
                      const long int* m,                                \
                      const long int* n,                                \
                      T* a,                                             \
                      const long int* lda,                              \
                      T* s,                                             \
                      const T* u,                                       \
                      const long int* ldu,                              \
                      const T* vt,                                      \
                      const long int* ldvt,                             \
                      long int* info,                                   \
                      nt2::details::workspace<T> & w)                   \
    {                                                                   \
      NT2_F77NAME( NAME )(jobu, jobvt, m, n, a, lda, s, u, ldu,         \
                          vt, ldvt, w.getw(), w.query(), info);         \
      w.resizew(w.neededsize());                                        \
      NT2_F77NAME( NAME )(jobu, jobvt, m, n, a, lda, s, u, ldu,         \
                          vt,ldvt,w.getw(),&w.neededsize(),info);       \
    }                                                                   \
    inline void gesvd(const char* jobu,                                 \
                      const char* jobvt,                                \
                      const long int* m,                                \
                      const long int* n,                                \
                      T* a,                                             \
                      const long int* lda,                              \
                      T* s,                                         \
                      const T* u,                                       \
                      const long int* ldu,                              \
                      const T* vt,                                      \
                      const long int* ldvt,                             \
                      long int* info)                                   \
    {                                                                   \
      nt2::details::workspace<T> w;                                     \
      gesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, info, w);   \
    }                                                                   \

    NT2_GESVD(sgesvd, float)
    NT2_GESVD(dgesvd, double)
#undef NT2_GESVD
     
  }
}


#endif

