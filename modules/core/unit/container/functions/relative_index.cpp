/*******************************************************************************
 *         Copyright 2003 & onward LASMEA UMR 6602 CNRS/Univ. Clermont II
 *         Copyright 2009 & onward LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#define NT2_UNIT_MODULE "nt2 container relative_index"

#include <nt2/table.hpp>
#include <nt2/include/functions/extent.hpp>
#include <nt2/include/functions/of_size.hpp>
#include <nt2/include/functions/relative_index.hpp>

#include <nt2/sdk/unit/module.hpp>
#include <nt2/sdk/unit/tests/basic.hpp>
#include <nt2/sdk/unit/tests/relation.hpp>

NT2_TEST_CASE( integral_subscript )
{
  using nt2::table;
  using nt2::of_size;
  typedef double T;
  boost::dispatch::meta::as_< boost::dispatch::meta::as_integer<T>::type > tgt;

  table<T> a( of_size(5,4,3,2) );

  for(int i=1;i<=5*4*3*2;i++)
    NT2_TEST_EQUAL( nt2::relative_index
                    ( boost::proto::child_c<1>(a(i,nt2::_)), 1, 5, 1, tgt )
                  , i
                  );
}

NT2_TEST_CASE( colon_subscript )
{
  using nt2::table;
  using nt2::of_size;
  typedef double T;
  boost::dispatch::meta::as_< boost::dispatch::meta::as_integer<T>::type > tgt;

  table<T> a( of_size(5,4,3,2) );

  for(int i=1;i<=5*4*3*2;i++)
    NT2_TEST_EQUAL( nt2::relative_index ( boost::proto::child_c<1>(a(nt2::_))
                                        , 1, 5, i, tgt
                                        )
                  , i
                  );
}

NT2_TEST_CASE( unity_colon_subscript )
{
  using nt2::table;
  using nt2::of_size;
  typedef double T;
  boost::dispatch::meta::as_< boost::dispatch::meta::as_integer<T>::type > tgt;

  table<T> a( of_size(5,4) );

  for(int i=2;i<=4;i++)
    NT2_TEST_EQUAL( nt2::relative_index ( boost::proto::child_c<1>(a(nt2::_(2,4)))
                                        , 1, 5, i-1, tgt
                                        )
                  , i
                  );
}

NT2_TEST_CASE( strided_colon_subscript )
{
  using nt2::table;
  using nt2::of_size;
  typedef double T;
  boost::dispatch::meta::as_< boost::dispatch::meta::as_integer<T>::type > tgt;

  table<T> a( of_size(15,4) );

  for(int i=1;i<=5;i++)
    NT2_TEST_EQUAL( nt2::relative_index ( boost::proto::child_c<1>(a(nt2::_(2,3,14)))
                                        , 1, 5, i, tgt
                                        )
                  , 2+(i-1)*3
                  );
}
