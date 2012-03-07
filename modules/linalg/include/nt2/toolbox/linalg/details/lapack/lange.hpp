/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_DETAILS_LAPACK_LANGE_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_DETAILS_LAPACK_LANGE_HPP_INCLUDED
/*
**  purpose
**  =======
**
**  xlange  returns the value of the one norm,  or the frobenius norm, or
**  the  infinity norm,  or the  element of  largest absolute value  of a
**  DATA TYPE matrix a.
**
**  description
**  ===========
**
**  xlange returns the value
**
**     xlange = ( max(abs(a(i,j))), norm = 'm' or 'm'
**              (
**              ( norm1(a),         norm = '1', 'o' or 'o'
**              (
**              ( normi(a),         norm = 'i' or 'i'
**              (
**              ( normf(a),         norm = 'f', 'f', 'e' or 'e'
**
**  where  norm1  denotes the  one norm of a matrix (maximum column sum),
**  normi  denotes the  infinity norm  of a matrix  (maximum row sum) and
**  normf  denotes the  frobenius norm of a matrix (square root of sum of
**  squares).  note that  max(abs(a(i,j)))  is not a  matrix norm.
**
**  arguments
**  =========
**
**  norm    (input) char
**          specifies the value to be returned in xlange as described
**          above.
**
**  m       (input) long int
**          the number of rows of the matrix a.  m >= 0.  when m = 0,
**          xlange is set to zero.
**
**  n       (input) long int
**          the number of columns of the matrix a.  n >= 0.  when n = 0,
**          xlange is set to zero.
**
**  a       (input) DATA TYPE array, dimension (lda,n)
**          the m by n matrix a.
**
**  lda     (input) long int
**          the leading dimension of the array a.  lda >= max(m,1).
**
**
** =====================================================================
*/
namespace nt2
{
  namespace details
  {
    //////////////////////////////////////////////////////////////////////
    // lange calls
    //////////////////////////////////////////////////////////////////////
    extern "C"
    {
      #define NT2_COMPLEX void
      float NT2_F77NAME(clange)(const char* norm, const long int* m, const long int* n,
                                const NT2_COMPLEX* a, const long int* lda, float* work);
      double NT2_F77NAME(dlange)(const char* norm, const long int* m, const long int* n,
                                 const double* a, const long int* lda, double* work);
      float NT2_F77NAME(slange)(const char* norm, const long int* m, const long int* n,
                                const float* a, const long int* lda, float* work);
      double NT2_F77NAME(zlange)(const char* norm, const long int* m, const long int* n,
                                 const NT2_COMPLEX* a, const long int* lda, double* work);
      #undef NT2_COMPLEX
    }

#define NT2_LANGE(NAME, T, TBASE)                               \
  inline TBASE lange(const char* norm,                          \
                     const long int* m,                         \
                     const long int* n,                         \
                     const T* a,                                \
                     const long int* lda,                       \
                     nt2::details::workspace<T> & w)            \
  {                                                             \
    w.resizerw((((*norm) == 'i')||((*norm) == 'I'))?*n:0);      \
    return NT2_F77NAME( NAME )(norm, m, n, a, lda, w.getrw());      \
  }                                                             \
  inline TBASE lange(const char* norm,                          \
                     const long int* m,                         \
                     const long int* n,                         \
                     const T* a,                                \
                     const long int* lda)                       \
  {                                                             \
    nt2::details::workspace<T> w;                               \
    return lange(norm, m, n, a, lda, w);                        \
  }                                                             \
    
  NT2_LANGE(slange, float, float)
  NT2_LANGE(dlange, double, double)
  NT2_LANGE(clange, std::complex<float> , float)
  NT2_LANGE(zlange, std::complex<double>, double)
    
#undef NT2_LANGE
    
    }  
}


#endif

// /////////////////////////////////////////////////////////////////////////////
// End of lange.hpp
// /////////////////////////////////////////////////////////////////////////////
