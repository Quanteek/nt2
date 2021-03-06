//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#ifndef NT2_CORE_FUNCTIONS_SCALAR_PARTIAL_FOLD_HPP_INCLUDED
#define NT2_CORE_FUNCTIONS_SCALAR_PARTIAL_FOLD_HPP_INCLUDED

#include <nt2/core/functions/partial_fold.hpp>
#include <boost/fusion/include/pop_front.hpp>
#include <nt2/include/functions/numel.hpp>
#include <nt2/sdk/details/type_id.hpp>


namespace nt2 { namespace details
{
  template <class X, class N, class B, class U>
  BOOST_FORCEINLINE typename X::value_type
  partial_fold_step(X const& in, const std::size_t& p, const std::size_t& mbound, N const& neutral, B const& bop, U const& uop)
  {
    typedef typename X::value_type   value_type;
    typedef typename X::extent_type  extent_type;
    extent_type ext = in.extent();

    std::size_t ibound  = boost::fusion::at_c<0>(ext);
    value_type out = neutral(nt2::meta::as_<value_type>());
      
    for(std::size_t m = 0, m_ = 0; m < mbound; ++m, m_+=ibound)
    {
        out = bop(out, nt2::run(in, m_+p, meta::as_<value_type>()));
    }

    
    return out;
  }
} }


namespace nt2 { namespace ext
{
  //============================================================================
  // Generates partial_fold
  //============================================================================
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::partial_fold_, tag::cpu_, (A0)(A1)(A2)(A3)(A4)
                            , (ast_< A0>)
                              (ast_< A1>)
                              (unspecified_<A2>)
                              (unspecified_<A3>)
                              (unspecified_<A4>)
                            )
  {
    typedef void                                                               result_type;
    typedef typename A0::value_type                                            value_type;
    typedef typename boost::remove_reference<A1>::type::extent_type            extent_type;

    BOOST_FORCEINLINE result_type operator()(A0& out, A1& in, A2 const& neutral, A3 const& bop, A4 const& uop) const
    {
      extent_type ext = in.extent();

      std::size_t ibound  = boost::fusion::at_c<0>(ext);
      std::size_t mbound =  boost::fusion::at_c<1>(ext);
      std::size_t obound =  boost::fusion::at_c<2>(ext);
      std::size_t id;

      for(std::size_t o = 0, o_ = 0; o < obound; ++o, o_+=ibound){
        for(std::size_t i = 0; i < ibound; ++i){
          id = i+o_;
          nt2::run(out, id, details::partial_fold_step(in,id,mbound,neutral,bop,uop));
        }
      }
    }

  };
    
  } }

#endif
