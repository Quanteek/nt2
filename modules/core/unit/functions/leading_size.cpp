/*******************************************************************************
 *         Copyright 2003 & onward LASMEA UMR 6602 CNRS/Univ. Clermont II
 *         Copyright 2009 & onward LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#define NT2_UNIT_MODULE "nt2::leading_size function"

#include <nt2/table.hpp>
#include <nt2/include/functions/leading_size.hpp>
#include <nt2/include/functions/of_size.hpp>

#include <nt2/sdk/unit/module.hpp>
#include <nt2/sdk/unit/tests/basic.hpp>
#include <nt2/sdk/unit/tests/relation.hpp>

////////////////////////////////////////////////////////////////////////////////
// leading_size of arithmetic types
////////////////////////////////////////////////////////////////////////////////
// NT2_TEST_CASE( fundamental_leading_size )
// {
//   using nt2::leading_size;

//   NT2_TEST_EQUAL( leading_size('4'), 1U  );
//   NT2_TEST_EQUAL( leading_size(4)  , 1U  );
//   NT2_TEST_EQUAL( leading_size(4.) , 1U  );
//   NT2_TEST_EQUAL( leading_size(4.f), 1U  );
// }

////////////////////////////////////////////////////////////////////////////////
// leading_size of container
////////////////////////////////////////////////////////////////////////////////
// NT2_TEST_CASE( container_leading_size )
// {
//   using nt2::leading_size;
//   using nt2::of_size;

//   typedef nt2::memory::container< nt2::tag::table_
//                                 , nt2::id_<0>,float,nt2::settings()
//                                 > container_t;

//   container_t t0;
//   container_t t1( of_size(2,1,1,1) );
//   container_t t2( of_size(2,2,1,1) );
//   container_t t3( of_size(2,2,2,1) );
//   container_t t4( of_size(2,2,2,2) );

//   NT2_TEST_EQUAL( size_t(leading_size(t0)), 0U   );
//   NT2_TEST_EQUAL( size_t(leading_size(t1)), 2U   );
//   NT2_TEST_EQUAL( size_t(leading_size(t2)), 4U   );
//   NT2_TEST_EQUAL( size_t(leading_size(t3)), 8U   );
//   NT2_TEST_EQUAL( size_t(leading_size(t4)), 16U  );
// }

////////////////////////////////////////////////////////////////////////////////
// leading_size of table
////////////////////////////////////////////////////////////////////////////////
NT2_TEST_CASE( table_leading_size )
{
  using nt2::leading_size;
  using nt2::of_size;
  using nt2::table;

  table<float> t0;
  table<float> t2( of_size(2) );
  table<float> t3( of_size(3) );
  table<float> t4( of_size(4) );
  table<float> t5( of_size(5) );
  table<float> t22( of_size(2,2) );
  table<float> t32( of_size(3,2) );
  table<float> t42( of_size(4,2) );

  NT2_TEST_EQUAL( size_t(leading_size(t0)), 0U   );
  NT2_TEST_EQUAL( size_t(leading_size(t2)), 2U   );
  NT2_TEST_EQUAL( size_t(leading_size(t3)), 3U   );
  NT2_TEST_EQUAL( size_t(leading_size(t4)), 4U  );
  NT2_TEST_EQUAL( size_t(leading_size(t5)), 5U  );  
  NT2_TEST_EQUAL( size_t(leading_size(t22)), 16U   );
  NT2_TEST_EQUAL( size_t(leading_size(t32)), 16U   );
  NT2_TEST_EQUAL( size_t(leading_size(t42)), 16U  );
}

////////////////////////////////////////////////////////////////////////////////
// leading_size of table expression
////////////////////////////////////////////////////////////////////////////////
// NT2_TEST_CASE( expression_leading_size )
// {
//   using nt2::leading_size;
//   using nt2::of_size;
//   using nt2::table;

//   table<float> t0;
//   table<float> t1( of_size(2) );
//   table<float> t2( of_size(2,2) );
//   table<float> t3( of_size(2,2,2) );
//   table<float> t4( of_size(2,2,2,2) );

//   NT2_TEST_EQUAL( size_t(leading_size(-t0)), 0U   );
//   NT2_TEST_EQUAL( size_t(leading_size(t1*t1)), 4U   );
//   NT2_TEST_EQUAL( size_t(leading_size(t2-t2*t2)), 4U   );
//   NT2_TEST_EQUAL( size_t(leading_size(t3/t3+t3)), 4U   );
//   NT2_TEST_EQUAL( size_t(leading_size(t4 * -t4)), 4U  );
// }
