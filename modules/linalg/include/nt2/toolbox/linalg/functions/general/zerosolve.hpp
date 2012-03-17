/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_GENERAL_SUBSPACE_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_GENERAL_SUBSPACE_HPP_INCLUDED
#include <nt2/include/functions/subspace.hpp>
#include <nt2/include/functions/svd.hpp>
#include <nt2/include/functions/abs.hpp>
#include <nt2/include/functions/asin.hpp>
#include <nt2/include/functions/min.hpp>
#include <nt2/include/functions/width.hpp>
#include <nt2/include/functions/height.hpp>
#include <nt2/include/functions/orth.hpp>
#include <nt2/include/functions/norm.hpp>
#include <nt2/include/constants/mone.hpp>
#include <nt2/table.hpp>

namespace nt2
{
    BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::subspace_, tag::cpu_,
                                       (A)(SA)(B), 
                                       ((expr_< table_<unspecified_<A>,SA>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                                       (scalar_<floating<B> > )
                                       )
  {
    typedef size_t result_type; 
    BOOST_DISPATCH_FORCE_INLINE result_type operator()(const A& a, const B epsi) const
    {
      return svd<A>(a).zerosolve(epsi); 
    }
  };
  
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::subspace_, tag::cpu_,
                                     (A)(SA), 
                                     ((expr_< table_<unspecified_<A>,SA>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                                     )
  {
    typedef typename A:value_type value_type;
    typedef A result_type; 
    BOOST_DISPATCH_FORCE_INLINE result_type operator()(const A& a) const
    {
      return svd<A>(a).zerosolve(); 
    }
  };

}


#endif

