//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#ifndef NT2_CORE_UTILITY_GENERATIVE_CATEGORY_HPP_INCLUDED
#define NT2_CORE_UTILITY_GENERATIVE_CATEGORY_HPP_INCLUDED

#include <nt2/sdk/meta/property_of.hpp>
#include <nt2/sdk/meta/hierarchy_of.hpp>

namespace boost { namespace dispatch { namespace meta
{
  template<class T, class F>
  struct generative_ : generative_< typename T::parent, F >
  {
    typedef generative_<typename T::parent, F> parent;
  };

  template<class T, class F>
  struct  generative_< unspecified_<T>, F >
        : generic_< typename property_of<T>::type >
  {
    typedef generic_< typename property_of<T>::type > parent;
  };
} } }

namespace nt2 { namespace ext
{
  using boost::dispatch::meta::generative_;
} }

#endif
