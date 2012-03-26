//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#ifndef NT2_CORE_FUNCTIONS_SCALAR_FOLD_HPP_INCLUDED
#define NT2_CORE_FUNCTIONS_SCALAR_FOLD_HPP_INCLUDED

#include <nt2/core/functions/fold.hpp>

namespace nt2 { namespace ext
{
  //============================================================================
  // Generates fold
  //============================================================================
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::fold_, tag::cpu_, (A1)(A2)(A3)
                              , ((ast_< A1>))
                              (unspecified_<A2>)
                              (unspecified_<A3>)
                            )
  {
    typedef typename boost::remove_reference<A1>::type::extent_type            extent_type;
    typedef typename boost::remove_reference<A1>::type::value_type             result_type;

    BOOST_FORCEINLINE result_type operator()(A1& in, A2 const& neutral, A3 const& op ) const
    {
      extent_type ext = in.extent();
      result_type out = neutral(nt2::meta::as_<result_type>());
      for(std::size_t c_0 = 1; c_0 <=ext[0]; ++c_0){
        out = op(out, in(c_0));
      }
      return out;     

    }
  };

} }

#endif
