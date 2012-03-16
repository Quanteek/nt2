/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_GENERAL_NORM_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_GENERAL_NORM_HPP_INCLUDED
#include <nt2/include/functions/norm.hpp>
#include <nt2/include/functions/isvector.hpp>
#include <nt2/include/functions/vecnorm.hpp>
#include <nt2/include/functions/is_nan.hpp>
#include <nt2/include/functions/is_finite.hpp>
#include <nt2/include/functions/is_gtz.hpp>
#include <nt2/include/functions/pow.hpp>
#include <nt2/include/functions/rec.hpp>
#include <nt2/include/functions/globalMax.hpp>
#include <nt2/toolbox/linalg/details/lapack/lange.hpp>
#include <nt2/include/functions/height.hpp>
#include <nt2/include/functions/width.hpp>
#include <nt2/include/functions/leading_size.hpp>

namespace nt2
{
    BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::norm_, tag::cpu_,
                                       (A0)(SA)(C), 
                                       ((expr_< table_<unspecified_<A0>,SA>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                                       (scalar_<aritmetic_<C> > )
                                       )
  {
    typedef typename A:value_type result_type; 
    BOOST_SIMD_FUNCTOR_CALL(2)
    {
      if (is_vector(a0))
        {
          return vecnorm(a0, a1); 
        }
      else if (is_matrix(a0))
        {
          la_int m = height(a0), n = width(a0);
          la_int lda0 = leading_size(a0); 
          if (a1 == 0){ 
            const char c = 'I';
            return nt2::details::lange(&c, &m, &n, (base_t*)const_cast<A0&>(a0).raw(), &lda0);
          } else if (p == 2){
            return svd<A0>(a0).get_w()(0); 
          } else if (p == 1) {
            char c = '1'; 
            return nt2::details::lange(&c, &m, &n, (base_t*)const_cast<A0&>(l1).raw(), &lda0);
          } else if (p == -1){
            char c = 'F'; 
            return nt2::details::lange(&c, &m, &n, (base_t*)const_cast<A0&>(l1).raw(), &lda0);
          } else if (!isfinite(a1)){
            return globalMax(sum(abs(trans(l1))));         
          } else {
            assert(false, "Sorry Not Yet Implemented"); 
          }
        }
      else
        {
          BOOST_ASSERT_MSG(false, "a0 is not matrix or vector"); 
          return 0; 
        }
    }
  };

  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::norm_, tag::cpu_,
                                       (A0)(SA)(C), 
                                       ((expr_< table_<unspecified_<A0>,SA>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                                       (unspecified_<C> )
                                       )
  {
    typedef typename A:value_type result_type; 
    BOOST_DISPATCH_FORCE_INLINE result_type operator()(A& a, const char * a1) const
    {
      if (strcmp(std::string(a1).c_str(), "frob") == 0)
        {
          return norm(a0, -1);
        }
      else if (strcmp(std::string(a1).c_str(), "inf") == 0)
        {
          return norm(a0, 0);
        }
      else
        {
          BOOST_ASSERT_MSG(false, "unknown option"); 
          return 0; 
        }
    }
  }
};

}


#endif

