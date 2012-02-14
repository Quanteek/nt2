//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II       
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI       
//                                                                              
//          Distributed under the Boost Software License, Version 1.0.          
//                 See accompanying file LICENSE.txt or copy at                 
//                     http://www.boost.org/LICENSE_1_0.txt                     
//==============================================================================
#ifndef NT2_TOOLBOX_ALGEBRA_FUNCTIONS_BLAS_GENERAL_GEMV_HPP_INCLUDED
#define NT2_TOOLBOX_ALGEBRA_FUNCTIONS_BLAS_GENERAL_GEMV_HPP_INCLUDED

#include <nt2/table.hpp>
#include <nt2/core/container/category.hpp>
#include <nt2/include/functions/size.hpp>
#include <nt2/toolbox/algebra/blas/blas2.hpp>
#include <boost/simd/sdk/memory/align_on.hpp>
#include <nt2/include/functions/gemm.hpp>

#define F77NAME(x) x##_

namespace nt2
{
  namespace ext
  {
    
#define NT2_GEMV(T, PREFIX)                                             \
    inline void gemv(const char *ta, const long int *m,                 \
                     const long int *n,                                 \
                     const T *al,                                       \
                     const T *a, const long int *lda,                   \
                     const T *x, const long int *incx,                  \
                     const T *be,                                       \
                     T *y, const long int *incy)                        \
    {                                                                   \
      F77NAME(PREFIX ## gemv)(ta,m,n,al,a,lda,x,incx,be,y, incy);       \
    }                                                                   \
        
    NT2_GEMV(double, d)
      NT2_GEMV(float,  s)
      NT2_GEMV(std::complex<double>, z)
      NT2_GEMV(std::complex<float>, c)
#undef NT2_GEMV
      
    /////////////////////////////////////////////////////////////////////////////
    // Implementation when type A0 is table_<double_>
    /////////////////////////////////////////////////////////////////////////////
      NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::general_gemv_, tag::cpu_, 
                                  (A0)(S0)(A1)(S1)(A2)(S2)(A3)(A4)(A5), //(A6)
                                  ((table_< floating_<A0>, S0 > ))
                                  ((table_< floating_<A1>, S1 > ))
                                  ((table_< floating_<A2>, S2 > ))
                                  (unspecified_ < A5 > ) //TO DO Specify
                                  (scalar_ < arithmetic_<A3 > > )
                                  (scalar_ < arithmetic_<A4 > > )
                                  // (scalar_ < integer_<A5 > > )
                                  // (scalar_ < integer_<A6 > > )
                                  )
    {
      typedef int result_type;
      BOOST_FORCEINLINE result_type operator()(A0& y, A1 const& a, A2 const& x,
                                               A5 const& transa, 
                                               A3 const& al, A4 const& be//, 
                                               //    A5 const& inx, A6 const& iny
                                               )
      {
        //       assert(is_vector(x));
        //       assert(is_vector(y));
        typedef typename A0::value_type value_type; 
        const long int m   = nt2::size(a)(transa=='T'?2:1);
        const long int n   = nt2::size(a)(transa=='T'?1:2);
        const value_type alpha = al; 
        const long int lda = nt2::details::padding(a);
        const value_type beta  = be;
        const long int incx = 1; //inx;
        const long int incy = 1; //iny;
        //        const char transa = 'N'; 
        gemv(&transa,&m,&n,&alpha,a.begin(),&lda,x.begin(),&incx,&beta,y.begin(),&incy);
        return 0; 
      }
    };
    
  }
  template < char T,  class A0,  class A1,  class A2,  class A3,  class A4/*,  class A5,  class A6*/>
  inline void gemv(A0& a0, A1 const& a1, A2 const& a2,A3 const& a3,A4 const&  a4/*,A5 const& a5 = 1,A6 const&  a6 = 1*/)
  {
    general_gemv(a0, a1, a2, 'N', a3, a4); //, a5, a6); 
  }
  
  template < char T,  class A0,  class A1,  class A2,  class A3>
  inline void gemv(A0& a0, A1 const& a1, A2 const& a2,A3 const& a3)
  {
    typedef typename A0::value_type value_type; 
    general_gemv(a0, a1, a2, 'N', a3, Zero<value_type>()); //, 1, 1); 
  }
  
  template < char T,  class A0,  class A1,  class A2>
  inline void gemv(A0& a0, A1 const& a1, A2 const& a2)
  {
    typedef typename A0::value_type value_type; 
    general_gemv(a0, a1, a2, 'N', One<value_type>(), Zero<value_type>()); //, 1, 1); 
  }
}

#undef F77NAME

#endif
