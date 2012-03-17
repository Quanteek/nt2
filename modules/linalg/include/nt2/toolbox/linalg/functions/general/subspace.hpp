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
#include <nt2/include/functions/abs.hpp>
#include <nt2/include/functions/asin.hpp>
#include <nt2/include/functions/min.hpp>
#include <nt2/include/functions/width.hpp>
#include <nt2/include/functions/height.hpp>
#include <nt2/include/functions/orth.hpp>
#include <nt2/include/functions/norm.hpp>
#include <nt2/table.hpp>

namespace nt2
{
    BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::subspace_, tag::cpu_,
                                       (A)(SA)(B)(SB), 
                                       ((expr_< table_<unspecified_<A>,SA>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                                       ((expr_< table_<unspecified_<B>,SB>,nt2::tag::terminal_,boost::mpl::long_<0> >))
                                       )
  {
    typedef typename A:value_type result_type;
    BOOST_DISPATCH_FORCE_INLINE result_type operator()(const A& a, const B& b) const
    {
      typedef table<value_type> tab_t;
      bool t =  nt2::width(a) >=  nt2::width(b); 
      tab_t oa(t ? a : b);
      tab_t ob(t ? b : a);
      oa = orth(oa);
      ob = orth(ob);
      for (size_t k=0; k < nt2::width(oa);  ++k){
         ob -= oa(_,k)*prodtMM(conj(oa(_,k)), ob);
       }
       // Make sure it's magnitude is less than 1.
      return nt2::asin(nt2::min(1,(nt2::norm(ob, 2))));
     }
  };

}


#endif

