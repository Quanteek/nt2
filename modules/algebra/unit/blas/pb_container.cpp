//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II       
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI       
//                                                                              
//          Distributed under the Boost Software License, Version 1.0.          
//                 See accompanying file LICENSE.txt or copy at                 
//                     http://www.boost.org/LICENSE_1_0.txt                     
//==============================================================================
#define NT2_UNIT_MODULE "nt2 algebra toolbox - pb_container"

//////////////////////////////////////////////////////////////////////////////
// unit test behavior of algebra components
//////////////////////////////////////////////////////////////////////////////

#include <nt2/table.hpp>
#include <nt2/include/functions/gemv.hpp>
#include <nt2/include/functions/size.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/make_vector.hpp>
#include <nt2/sdk/unit/tests.hpp>
#include <nt2/sdk/unit/module.hpp>

NT2_TEST_CASE_TPL ( pb_container, (float)) 
{
  using nt2::meta::make_container; 
  using nt2::of_size_;

  {
    //   // ok //
    //   typedef typename make_container<nt2::tag::table_, T, of_size_<17,5>  >::type table_type_a;
    //   table_type_a a;
    //   boost::fusion::vector<int,int> pos1;
    //   pos1 = boost::fusion::make_vector(1, 1);
    //   a[pos1] =22;
  } 
  {
    //   // ok //
    //   typedef typename make_container<nt2::tag::table_, T, of_size_<5> >::type table_type_b;
    //   table_type_b b;
    //   boost::fusion::vector<int> pos2;
    //   pos2 = boost::fusion::make_vector(1);
    //   b[pos2] =22;
  } 
  {
    //   // ok //
    //   typedef typename make_container<nt2::tag::table_, T, of_size_<1, 5> >::type table_type_c;
    //   table_type_c c;
    //   boost::fusion::vector<int, int> pos3;
    //   pos3 = boost::fusion::make_vector(1, 1);
    //   c[pos3] =22;
  } 
  {
    // not ok //
    typedef typename make_container<nt2::tag::table_, T, of_size_<5, 1> >::type table_type_c;
    table_type_c c;
    boost::fusion::vector<int, int> pos3;
    pos3 = boost::fusion::make_vector(1, 1);
    c[pos3] =22;
  } 
  {
    // not ok //
    typedef typename make_container<nt2::tag::table_, T, of_size_<5> >::type table_type_d;
    table_type_d d;
    boost::fusion::vector<int, int> pos4;
    pos4 = boost::fusion::make_vector(1, 1);
    d[pos4] =22;
  } 
}
