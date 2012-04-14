/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_GENERAL_VECNORM_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_GENERAL_VECNORM_HPP_INCLUDED
#include <nt2/include/functions/vecnorm.hpp>
#include <nt2/include/functions/isvector.hpp>
#include <nt2/include/functions/norm_eucl.hpp>
#include <nt2/include/functions/globalNormp.hpp>
#include <nt2/include/functions/is_nan.hpp>
#include <nt2/include/functions/is_finite.hpp>
#include <nt2/include/functions/is_gtz.hpp>
#include <nt2/include/functions/pow.hpp>
#include <nt2/include/functions/rec.hpp>

namespace nt2
{
    BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::vecnorm_, tag::cpu_,
                                       (A)(SA)(C), 
                                       ((expr_< table_<unspecified_<A>,SA>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                                       (scalar_<arithmetic_<C> > )
                                       )
  {
    typedef typename A:value_type result_type; 
    BOOST_SIMD_FUNCTOR_CALL(2)
    {
      if (is_nan(a1) return Nan<value_type>(); 
      if (a1 == 2){
        return norm_eucl(a0(_));
      } else if (nt2::is_finite(a1)){
        return nt2::pow(normp(a1(_), p), nt2::rec(result_type(p)));
      } else if (is_gtz(a1)){
        return nt2::max(abs(a0(_)));
      } else {
        return nt2::min(abs(a0(_)));
      }
    }
  };

}


#endif

