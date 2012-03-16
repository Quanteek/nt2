/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_GENERAL_NORMEST_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_GENERAL_NORMEST_HPP_INCLUDED
#include <nt2/include/functions/normest.hpp>
#include <nt2/include/functions/abs.hpp>
#include <nt2/include/functions/is_eqz.hpp>
#include <nt2/include/constants/zero.hpp>
#include <nt2/include/constants/sqrteps.hpp>
// #include <nt2/include/functions/isvector.hpp>
// #include <nt2/include/functions/vecnormest.hpp>
// #include <nt2/include/functions/is_nan.hpp>
// #include <nt2/include/functions/is_finite.hpp>
// #include <nt2/include/functions/is_gtz.hpp>
// #include <nt2/include/functions/pow.hpp>
// #include <nt2/include/functions/rec.hpp>
// #include <nt2/include/functions/globalMax.hpp>
// #include <nt2/toolbox/linalg/details/lapack/lange.hpp>
// #include <nt2/include/functions/height.hpp>
// #include <nt2/include/functions/width.hpp>
// #include <nt2/include/functions/leading_size.hpp>

namespace nt2
{
    BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::normest_, tag::cpu_,
                                       (A0)(SA)(C), 
                                       ((expr_< table_<unspecified_<A0>,SA>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                                       (scalar_<floating_<C> > )
                                       )
  {
    typedef typename A:value_type result_type;
    BOOST_DISPATCH_FORCE_INLINE result_type operator()(const A0& s, const result_type &tol) const
    {
      typedef table<value_type> tab_t; 
      tab_t x =  nt2::sum(nt2::abs(s));
      size_t cnt = 0;
      result_type e =  globalNorm_eucl(x);
      if (is_eqz(e)) return e;
      x /=  e;
      result_type e0 =  Zero<result_type>();
       while (nt2::abs(e-e0) > tol*e)
         {
           e0 = e;
           tab_t  sx =  x*trans(s);
           e =  globalNorm_eucl(sx);
           x = sx*s; ;
           if (all(nt2::abs(sx))== 0) sx = rand(size(sx), meta::as<result_type>());
           x /= globalNorm_eucl(x);
           if (++cnt > 100) break; 
         }
       return e;     
    }
  };

  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::normest_, tag::cpu_,
                                       (A0)(SA)(C), 
                                       ((expr_< table_<unspecified_<A0>,SA>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                                       )
  {
    typedef typename A:value_type result_type; 
    BOOST_DISPATCH_FORCE_INLINE result_type operator()(A& a) const
    {
      return nt2::normest(a, Sqrteps<result_type>()); 
    }
  }
};

}


#endif

