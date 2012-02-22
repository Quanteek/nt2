//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II       
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI       
//                                                                              
//          Distributed under the Boost Software License, Version 1.0.          
//                 See accompanying file LICENSE.txt or copy at                 
//                     http://www.boost.org/LICENSE_1_0.txt                     
//==============================================================================
#ifndef NT2_TOOLBOX_ALGEBRA_FUNCTIONS_BLAS_TRMV_HPP_INCLUDED
#define NT2_TOOLBOX_ALGEBRA_FUNCTIONS_BLAS_TRMV_HPP_INCLUDED

#include <nt2/table.hpp>
#include <nt2/core/container/category.hpp>
#include <nt2/include/functions/size.hpp>
#include <nt2/include/functions/extent.hpp>
#include <nt2/include/functions/numel.hpp>
#include <nt2/toolbox/algebra/blas/blas2.hpp>
#include <nt2/toolbox/algebra/details/padding.hpp>
#include <nt2/include/constants/zero.hpp>
#include <nt2/include/constants/one.hpp>
#include <nt2/include/functions/isvector.hpp>
#include <nt2/include/functions/issquare.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/assert.hpp>

namespace nt2
{
  namespace ext
  {
    
#define NT2_TRMV(T, PREFIX)                                             \
    inline void trmv(const char *uplo, const char *trans,               \
                     const char *diag,                                  \
                     const long int *n,                                 \
                     const T *a,                                        \
                     const long int *lda,                               \
                     T *x, const long int *incx)                        \
    {                                                                   \
      BOOST_PP_CAT(PREFIX,BOOST_PP_CAT(trmv,_))(uplo,trans,diag,n,a,lda,x,incx);       \
    }                                                                   \
        
    NT2_TRMV(double, d)
    NT2_TRMV(float,  s)
    NT2_TRMV(std::complex<double>, z)
    NT2_TRMV(std::complex<float>, c)

#undef NT2_TRMV
      
    /////////////////////////////////////////////////////////////////////////////
    // Implementation when table_ are floating_
    /////////////////////////////////////////////////////////////////////////////
      NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::trmv_, tag::cpu_
                                , (A0)(S0)(A1)(S1)(ST)
                                , (unspecified_ < ST > )
                                  ((table_< floating_<A0>, S0 > ))
                                  ((table_< floating_<A1>, S1 > ))
                                  )
    {
      typedef void result_type;
      BOOST_FORCEINLINE result_type operator()( ST const&, 
                                                A0 const& a,
                                                A1 & x
                                                )
      {
        const char upload = ST::upload;
        const char transa = ST::transa;
        const char transx = ST::transx;
        const char diag   = ST::diag; 
        typedef typename A0::value_type value_type; 

        //These are not to be asserted
        //         BOOST_ASSERT_MSG( ((diag == 'N')||nt2::is_unit(a)), 
        //                           "Matrix-vector product y = A°*x° (trmv) must be called with unit-triangular table A");
        //         BOOST_ASSERT_MSG( (nt2::is_triangular(a, upload)), 
        //                           "Matrix-vector product y = A°*x° (trmv) must be called with triangular table A");
        BOOST_ASSERT_MSG( (nt2::is_square(a)), 
                          "Matrix-vector product y = A°*x° (trmv) must be called with square table A");

        const long int n   = nt2::extent(a)[0]; 

        BOOST_ASSERT_MSG( (n == nt2::numel(x)), //TODO replace numel by proper dim here
                          "In matrix-vector product y = A°*x° (trmv) the inner number of columns and lines \
                           of A must be equal to the number of columns of x° ");
        
        const long int lda = nt2::details::padding(a);
        const long int incx = 1; //(transx == 'N') ? 1 : nt2::details::padding(x)/sizeof(value_type); 
        trmv(&upload,&transa,&diag,&n,a.begin(),&lda,x.begin(),&incx);
      }
    };
    
  }

}

#endif
