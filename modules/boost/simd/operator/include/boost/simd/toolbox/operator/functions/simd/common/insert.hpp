//==============================================================================
//         Copyright 2003 - 2011 LASMEA UMR 6602 CNRS/Univ. Clermont II         
//         Copyright 2009 - 2011 LRI    UMR 8623 CNRS/Univ Paris Sud XI         
//                                                                              
//          Distributed under the Boost Software License, Version 1.0.          
//                 See accompanying file LICENSE.txt or copy at                 
//                     http://www.boost.org/LICENSE_1_0.txt                     
//==============================================================================
#ifndef BOOST_SIMD_TOOLBOX_OPERATOR_FUNCTIONS_SIMD_COMMON_INSERT_HPP_INCLUDED
#define BOOST_SIMD_TOOLBOX_OPERATOR_FUNCTIONS_SIMD_COMMON_INSERT_HPP_INCLUDED

#include <boost/simd/toolbox/operator/functions/insert.hpp>
#include <boost/simd/include/constants/allbits.hpp>
#include <boost/simd/sdk/meta/scalar_of.hpp>
#include <boost/simd/sdk/details/aliasing.hpp>
#include <boost/simd/sdk/simd/logical.hpp>
#include <boost/mpl/equal_to.hpp>
#include <boost/mpl/sizeof.hpp>

#include <boost/mpl/assert.hpp>
#include <boost/fusion/include/at.hpp>
#include <boost/fusion/include/size.hpp>
#include <boost/simd/sdk/meta/iterate.hpp>

// workaround for circular includes
namespace boost { namespace simd { namespace tag
{
  struct Allbits;
} } }

namespace boost { namespace simd { namespace ext
{
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( boost::simd::tag::insert_, tag::cpu_, (A0)(A1)(A2)(X)
                            , (fusion_sequence_<A0>)
                              ((simd_< fusion_sequence_<A1>, X >))
                              (scalar_< integer_<A2> >)
                            )
  {
    struct insert_fusion
    {
      insert_fusion(A0 const& a0_, A1& a1_, A2 const& a2_)
        : a0(a0_), a1(a1_), a2(a2_)
      {
      }

      template<int I>
      void operator()() const
      {
        insert(fusion::at_c<I>(a0), fusion::at_c<I>(a1), a2);
      }

      A0 const& a0;
      A1& a1;
      A2 const& a2;
    };

    typedef A1& result_type;
    BOOST_FORCEINLINE result_type operator()(A0 const& a0, A1& a1, A2 const& a2) const
    {
      BOOST_MPL_ASSERT_MSG( fusion::result_of::size<A0>::type::value == fusion::result_of::size<A1>::type::value
                          , BOOST_SIMD_INSERT_FUSION_SEQUENCE_SIZE_MISMATCH
                          , (A0, A1)
                          );

      static const int N = fusion::result_of::size<A0>::type::value;
      meta::iterate<N>(insert_fusion(a0, a1, a2));
      return a1;
    }
  };

  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( boost::simd::tag::insert_, tag::cpu_, (A0)(A1)(A2)(X)
                            , (scalar_< arithmetic_<A0> >)
                              ((simd_< arithmetic_<A1>, X >))
                              (scalar_< integer_<A2> >)
                            )
  {
    typedef A1& result_type;
    BOOST_FORCEINLINE result_type operator()(A0 const& a0, A1& a1, A2 const& a2) const
    {
      typedef typename meta::scalar_of<A1>::type stype;
      reinterpret_cast<stype BOOST_SIMD_MAY_ALIAS*>(&a1)[a2] = a0;
      return a1;
    }
  };
  
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION_IF( boost::simd::tag::insert_, tag::cpu_, (A0)(A1)(A2)(X)
                            , (mpl::equal_to< mpl::sizeof_<A1>, mpl::sizeof_<typename A1::type> >)
                            , (scalar_< logical_<A0> >)
                              ((simd_< logical_<A1>, X >))
                              (scalar_< integer_<A2> >)
                            )
  {
    typedef A1& result_type;
    BOOST_FORCEINLINE result_type operator()(A0 const& a0, A1& a1, A2 const& a2) const
    {
      typedef typename meta::scalar_of<typename A1::type>::type type;
      type allbits = typename dispatch::make_functor<tag::Allbits, A0>::type()(dispatch::meta::as_<type>());
      insert(a0 ? allbits : a0, reinterpret_cast<typename A1::type&>(a1), a2);
      return a1;
    }
  };
  
} } }

#endif
