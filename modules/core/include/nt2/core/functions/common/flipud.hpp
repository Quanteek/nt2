//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#ifndef NT2_CORE_FUNCTIONS_COMMON_FLIPUD_HPP_INCLUDED
#define NT2_CORE_FUNCTIONS_COMMON_FLIPUD_HPP_INCLUDED

#include <nt2/core/functions/flipud.hpp>
#include <nt2/include/functions/ge.hpp>
#include <nt2/include/functions/run.hpp>
#include <nt2/include/constants/zero.hpp>
#include <nt2/include/functions/splat.hpp>
#include <nt2/include/functions/sub2ind.hpp>
#include <nt2/include/functions/ind2sub.hpp>
#include <nt2/include/functions/height.hpp>
#include <nt2/include/functions/if_else.hpp>
#include <nt2/include/functions/enumerate.hpp>

namespace nt2 { namespace ext
{
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::run_, tag::cpu_
                            , (A0)(State)(Data)(N)
                            , ((node_<A0, nt2::tag::flipud_, N>))
                              (generic_< integer_<State> >)
                              ((unspecified_<Data>))
                            )
  {
    typedef typename meta::strip<Data>::type::type                    result_type;
    typedef typename meta::as_integer<result_type>::type              i_t;
    typedef typename meta::call<nt2::tag::ind2sub_(_2D,State)>::type  sub_t;

    BOOST_FORCEINLINE result_type
    operator()(A0 const& a0, State const& p, Data const& t) const
    {
      // Retrieve 2D position from the linear index
      sub_t pos = ind2sub(_2D(a0.extent()),p);
      std::size_t n = height(a0); 
      return  nt2::run(boost::proto::child_c<0>(a0),p,t); //INCORRECT MUST CHANGE run here
    }
  };

} }

#endif
