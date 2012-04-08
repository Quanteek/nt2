/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_FUNCTIONS_CIRCUL_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_FUNCTIONS_CIRCUL_HPP_INCLUDED
#include <nt2/include/functions/circul.hpp>
#include <nt2/include/functions/colvect.hpp>
#include <nt2/include/functions/rowvect.hpp>
#include <nt2/include/functions/flipud.hpp>
#include <nt2/include/functions/toeplitz.hpp>
#include <nt2/table.hpp>

namespace nt2
{
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::circul_, tag::cpu_,
                                     (A)(SA), 
                                     ((expr_< table_<unspecified_<A0>,SA>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                                     )
  {
    typedef typename A:value_type value_type;
    typedef typename A:index_type index_type;
    typedef table<value_type, index_type> result_type; 
    inline result_type operator()(const A& v)const
      {
        return toeplitz(flipud(colvect(v)),rowvect(v));
      }
  };
}


#endif
