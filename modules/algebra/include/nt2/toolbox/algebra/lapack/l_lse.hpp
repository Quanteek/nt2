/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2009-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_ALGEBRA_LAPACK_L_LSE_HPP_INCLUDED
#define NT2_TOOLBOX_ALGEBRA_LAPACK_L_LSE_HPP_INCLUDED

#include <nt2/toolbox/algebra/blas/f77_wrapper.hpp>



namespace nt2
{
  namespace details
  {
#define NT2_GGLSE(NAME, T)						\
    inline void gglse( const long int* m,				\
		       const long int* n,				\
		       const long int* p,				\
		       T* a, const long int* lda,			\
		       T* b, const long int* ldb,			\
		       T* c,						\
		       const T* d,					\
		       T* x,						\
		       long int* info,					\
		       workspace<T> & w)				\
    {									\
      NT2_F77NAME(NAME)(m, n, p, a, lda, b, ldb,			\
			c, d, x, w.getw(), w.query(), info);		\
      w.resizew(w.neededsize());					\
      NT2_F77NAME(NAME)(m, n, p, a, lda, b, ldb,			\
			c, d, x, w.getw(), &w.neededsize(), info);	\
    }									\
    inline void gglse( const long int* m,				\
		       const long int* n,				\
		       const long int* p,				\
		       T* a, const long int* lda,			\
		       T* b, const long int* ldb,			\
		       T* c,						\
		       const T* d,					\
		       T* x,						\
		       long int* info)					\
    {									\
      workspace<T> w;							\
      gglse(m, n, p, a, lda, b, ldb, c, d, x, info, w);			\
    }									\
	
    NT2_GGLSE(sgglse, float)
    NT2_GGLSE(dgglse, double)
    NT2_GGLSE(sgglse, std::complex<float>)
    NT2_GGLSE(dgglse, std::complex<double>)

#undef NT2_GGLSE
    
  }
}


#endif

// /////////////////////////////////////////////////////////////////////////////
// End of l_lse.hpp<3>
// /////////////////////////////////////////////////////////////////////////////
