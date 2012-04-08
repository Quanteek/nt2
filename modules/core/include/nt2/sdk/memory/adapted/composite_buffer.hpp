//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#ifndef NT2_SDK_META_MEMORY_ADAPTED_ARRAY_BUFFER_HPP_INCLUDED
#define NT2_SDK_META_MEMORY_ADAPTED_ARRAY_BUFFER_HPP_INCLUDED

#include <boost/mpl/apply.hpp>
#include <boost/dispatch/meta/model_of.hpp>
#include <boost/dispatch/meta/value_of.hpp>

namespace nt2 {  namespace memory
{
  //============================================================================
  // Forward declaration
  //============================================================================
  template<typename Buffer> struct composite_buffer;
} }

namespace boost { namespace dispatch { namespace meta
{
  //============================================================================
  // value_of specialization
  //============================================================================
  template<typename Buffer>
  struct value_of< nt2::memory::composite_buffer<Buffer> > : value_of<Buffer>
  {};

  //============================================================================
  // model_of specialization
  //============================================================================
  template<typename Buffer>
  struct model_of< nt2::memory::composite_buffer<Buffer> >
  {
    struct type
    {
      template<class X> struct apply
      {
        typedef typename  boost::mpl::
                          apply<typename model_of<Buffer>::type,X>::type  base_t;
        typedef nt2::memory::composite_buffer<base_t>                     type;
      };
    };
  };
} } }

#endif
