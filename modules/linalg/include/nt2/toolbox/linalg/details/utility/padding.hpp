//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II       
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI       
//                                                                              
//          Distributed under the Boost Software License, Version 1.0.          
//                 See accompanying file LICENSE.txt or copy at                 
//                     http://www.boost.org/LICENSE_1_0.txt                     
//==============================================================================
#ifndef NT2_TOOLBOX_LINALG_DETAILS_UTILITY_DETAILS_PADDING_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_DETAILS_UTILITY_DETAILS_PADDING_HPP_INCLUDED
//#include <boost/simd/sdk/memory/align_on.hpp>

namespace nt2 { namespace details 
{
  template < class T > long int padding(const T & a)
  {
//     typedef typename T::parent::lead_t lead_t_a;
//     long int tmp = boost::simd::memory::align_on(size(a, 1), lead_t_a::value);
//     std::cout << tmp << "  " << sizeof(typename T::value_type) << std::endl; ;
    return a.leading_size(); 
  }

} }

#endif
