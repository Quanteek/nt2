/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_GENERAL_ROT90_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_GENERAL_ROT90_HPP_INCLUDED
#include <nt2/include/functions/rot90.hpp>
#include <nt2/include/functions/trans.hpp>
#include <nt2/include/functions/flipud.hpp>
#include <nt2/include/functions/fliplr.hpp>
#include <nt2/table.hpp>

namespace nt2
{
    BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::rot90_, tag::cpu_,
                                       (A)(SA)(C), 
                                       ((expr_< table_<unspecified_<A>,SA>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                                       (scalar_<integer_<C> > )
                                       )
  {
    typedef typename A:value_type value_type;
    typedef table<value_type> result_type; 
    BOOST_SIMD_FUNCTOR_CALL(2)
    {
      k = ((k%4)+4)%4; 
      if(k == 1)  return flipud(trans(a0));
      if(k == 2)  return fliplr(flipud(a0));
      if(k == 3)  return trans(flipud(a0));
      return a0;
    }
  };

}


#endif

