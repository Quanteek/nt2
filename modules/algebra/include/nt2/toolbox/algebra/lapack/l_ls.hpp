/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2009-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_ALGEBRA_LAPACK_L_LS_HPP_INCLUDED
#define NT2_TOOLBOX_ALGEBRA_LAPACK_L_LS_HPP_INCLUDED



namespace nt2
{
  namespace details
  {
#define NT2_GELSY(NAME, T)			\
    inline void gelsy(const long int* m,	\
		      const long int* n,	\
		      const long int* nrhs,	\
		      T* a,			\
		      const long int* lda,	\
		      T* b,			\
		      const long int* ldb,	\
		      long int* jpvt,		\
		      const T* rcond,		\
		      long int* rank,		\
		      long int* info,		\
		      workspace<T> & w)					\
    {									\
      NT2_F77NAME( NAME )(m, n, nrhs, a, lda, b, ldb, jpvt,		\
			  rcond, rank, w.getw(), w.query(), info);	\
      w.resizew(w.neededsize());					\
      NT2_F77NAME( NAME )(m, n, nrhs, a, lda, b, ldb, jpvt,		\
			  rcond, rank, w.getw(), &w.neededsize(), info); \
    }									\
	inline void gelsy(const long int* m,				\
			  const long int* n,				\
			  const long int* nrhs,				\
			  T* a,						\
			  const long int* lda,				\
			  T* b,						\
			  const long int* ldb,				\
			  long int* jpvt,				\
			  const T* rcond,				\
			  long int* rank,				\
			  long int* info)				\
	{								\
	  workspace<T> w;						\
	  gelsy(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, info, w); \
	}								\
	
    NT2_GELSY(sgelsy, float)
    NT2_GELSY(dgelsy, double)
    NT2_GELSY(cgelsy, std::complex<float>,  float)
    NT2_GELSY(zgelsy, std::complex<double>, double)

#undef NT2_GELSY
  }
}


#endif

// /////////////////////////////////////////////////////////////////////////////
// End of l_ls.hpp
// /////////////////////////////////////////////////////////////////////////////
