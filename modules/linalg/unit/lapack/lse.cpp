//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II       
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI       
//                                                                              
//          Distributed under the Boost Software License, Version 1.0.          
//                 See accompanying file LICENSE.txt or copy at                 
//                     http://www.boost.org/LICENSE_1_0.txt                     
//==============================================================================
#define NT2_UNIT_MODULE "nt2 linalg toolbox - lse"

//////////////////////////////////////////////////////////////////////////////
// unit test behavior of algebra components
//////////////////////////////////////////////////////////////////////////////

#include <nt2/table.hpp>
//#include <nt2/include/functions/lse.hpp>
#include <nt2/include/functions/size.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/make_vector.hpp>
#include <nt2/sdk/unit/tests.hpp>
#include <nt2/sdk/unit/module.hpp>

NT2_TEST_CASE_TPL ( lse, NT2_REAL_TYPES) 
{
  using nt2::meta::make_container; 
  using nt2::of_size_;
   const int M = 5;
  typedef typename make_container<nt2::tag::table_, T, of_size_<M, 1>  >::type table_type_b;
  table_type_b b;
  std::size_t dim1_b = size(b)(1);
}

