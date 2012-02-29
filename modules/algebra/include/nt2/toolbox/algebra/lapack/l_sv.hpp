/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2009-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_ALGEBRA_LAPACK_L_SV_HPP_INCLUDED
#define NT2_TOOLBOX_ALGEBRA_LAPACK_L_SV_HPP_INCLUDED



namespace nt2
{
  namespace details
  {
#define NT2_GESV(NAME, T)					\
    inline void gesv(const long int* n,				\
		     const long int* nrhs,			\
		     T* a,					\
		     const long int* lda,			\
		     long int* ipiv,				\
		     T* b,					\
		     const long int* ldb,			\
		     long int* info)				\
    {								\
      NT2_F77NAME( NAME )(n, nrhs, a, lda, ipiv, b, ldb, info);	\
    }								\

    NT2_GESV(sgesv, float)
    NT2_GESV(dgesv, double)
    NT2_GESV(cgesv, std::complex<float>,  float)
    NT2_GESV(zgesv, std::complex<double>, double)

#undef NT2_GESV

#define NT2_POSV(NAME, T)			\
      inline void posv(const char* uplo,	\
		       const long int* n,	\
		       const long int* nrhs,	\
		       T* a,			\
		       const long int* lda,	\
		       T* b,			\
		       const long int* ldb,	\
		       long int* info)		\
      {								\
	F77NAME( NAME )(uplo, n, nrhs, a, lda, b, ldb, info);	\
      }								\
	  
    NT2_POSV(sposv, float)
    NT2_POSV(dposv, double)
    NT2_POSV(cposv, std::complex<float>,  float)
    NT2_POSV(zposv, std::complex<double>, double)

#undef NT2_POSV


  }
}


#endif

// /////////////////////////////////////////////////////////////////////////////
// End of l_sv.hpp
// /////////////////////////////////////////////////////////////////////////////
