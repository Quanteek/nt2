//////////////////////////////////////////////////////////////////////////////
///   Copyright 2003 and onward LASMEA UMR 6602 CNRS/U.B.P Clermont-Ferrand
///   Copyright 2009 and onward LRI    UMR 8623 CNRS/Univ Paris Sud XI
///
///          Distributed under the Boost Software License, Version 1.0
///                 See accompanying file LICENSE.txt or copy at
///                     http://www.boost.org/LICENSE_1_0.txt
//////////////////////////////////////////////////////////////////////////////
#define NT2_UNIT_MODULE "nt2 arithmetic toolbox - fma/scalar Mode"

//////////////////////////////////////////////////////////////////////////////
// unit test behavior of arithmetic components in scalar mode
//////////////////////////////////////////////////////////////////////////////
/// created by jt the 01/12/2010
/// 
#include <nt2/toolbox/arithmetic/include/functions/fma.hpp>
#include <nt2/include/functions/ulpdist.hpp>
#include <boost/type_traits/is_same.hpp>
#include <nt2/sdk/functor/meta/call.hpp>
#include <nt2/sdk/meta/as_integer.hpp>
#include <nt2/sdk/meta/as_floating.hpp>
#include <nt2/sdk/meta/as_signed.hpp>
#include <nt2/sdk/meta/upgrade.hpp>
#include <nt2/sdk/meta/downgrade.hpp>
#include <nt2/sdk/meta/scalar_of.hpp>
#include <boost/dispatch/meta/as_floating.hpp>
#include <boost/type_traits/common_type.hpp>
#include <nt2/sdk/unit/tests.hpp>
#include <nt2/sdk/unit/module.hpp>

#include <nt2/toolbox/constant/constant.hpp>


NT2_TEST_CASE_TPL ( fma_real__3_0,  NT2_REAL_TYPES)
{
  
  using nt2::fma;
  using nt2::tag::fma_;
  typedef typename nt2::meta::as_integer<T>::type iT;
  typedef typename nt2::meta::call<fma_(T,T,T)>::type r_t;
  typedef typename nt2::meta::scalar_of<r_t>::type ssr_t;
  typedef typename nt2::meta::upgrade<T>::type u_t;
  typedef typename boost::common_type<T,T,T>::type wished_r_t;


  // return type conformity test 
  NT2_TEST( (boost::is_same < r_t, wished_r_t >::value) );
  std::cout << std::endl; 
  double ulpd;
  ulpd=0.0;


  // specific values tests
  NT2_TEST_ULP_EQUAL(fma(nt2::Inf<T>(), nt2::Inf<T>(), nt2::Inf<T>()), nt2::Inf<T>(), 0);
  NT2_TEST_ULP_EQUAL(fma(nt2::Minf<T>(), nt2::Minf<T>(), nt2::Minf<T>()), nt2::Nan<T>(), 0);
  NT2_TEST_ULP_EQUAL(fma(nt2::Mone<T>(), nt2::Mone<T>(), nt2::Mone<T>()), nt2::Zero<T>(), 0);
  NT2_TEST_ULP_EQUAL(fma(nt2::Nan<T>(), nt2::Nan<T>(), nt2::Nan<T>()), nt2::Nan<T>(), 0);
  NT2_TEST_ULP_EQUAL(fma(nt2::One<T>(), nt2::One<T>(), nt2::One<T>()), nt2::Two<T>(), 0);
  NT2_TEST_ULP_EQUAL(fma(nt2::Two<T>(),nt2::Two<T>(),nt2::Three<T>()), nt2::Seven<T>(), 0);
  NT2_TEST_ULP_EQUAL(fma(nt2::Zero<T>(), nt2::Zero<T>(), nt2::Zero<T>()), nt2::Zero<T>(), 0);
} // end of test for floating_

NT2_TEST_CASE_TPL ( fma_unsigned_int__3_0,  NT2_UNSIGNED_TYPES)
{
  
  using nt2::fma;
  using nt2::tag::fma_;
  typedef typename nt2::meta::as_integer<T>::type iT;
  typedef typename nt2::meta::call<fma_(T,T,T)>::type r_t;
  typedef typename nt2::meta::scalar_of<r_t>::type ssr_t;
  typedef typename nt2::meta::upgrade<T>::type u_t;
  typedef typename boost::common_type<T,T,T>::type wished_r_t;


  // return type conformity test 
  NT2_TEST( (boost::is_same < r_t, wished_r_t >::value) );
  std::cout << std::endl; 
  double ulpd;
  ulpd=0.0;


  // specific values tests
  NT2_TEST_ULP_EQUAL(fma(nt2::One<T>(), nt2::One<T>(), nt2::One<T>()), nt2::Two<T>(), 0);
  NT2_TEST_ULP_EQUAL(fma(nt2::Two<T>(),nt2::Two<T>(),nt2::Three<T>()), nt2::Seven<T>(), 0);
  NT2_TEST_ULP_EQUAL(fma(nt2::Zero<T>(), nt2::Zero<T>(), nt2::Zero<T>()), nt2::Zero<T>(), 0);
} // end of test for unsigned_int_

NT2_TEST_CASE_TPL ( fma_signed_int__3_0,  NT2_INTEGRAL_SIGNED_TYPES)
{
  
  using nt2::fma;
  using nt2::tag::fma_;
  typedef typename nt2::meta::as_integer<T>::type iT;
  typedef typename nt2::meta::call<fma_(T,T,T)>::type r_t;
  typedef typename nt2::meta::scalar_of<r_t>::type ssr_t;
  typedef typename nt2::meta::upgrade<T>::type u_t;
  typedef typename boost::common_type<T,T,T>::type wished_r_t;


  // return type conformity test 
  NT2_TEST( (boost::is_same < r_t, wished_r_t >::value) );
  std::cout << std::endl; 
  double ulpd;
  ulpd=0.0;


  // specific values tests
  NT2_TEST_ULP_EQUAL(fma(nt2::Mone<T>(), nt2::Mone<T>(), nt2::Mone<T>()), nt2::Zero<T>(), 0);
  NT2_TEST_ULP_EQUAL(fma(nt2::One<T>(), nt2::One<T>(), nt2::One<T>()), nt2::Two<T>(), 0);
  NT2_TEST_ULP_EQUAL(fma(nt2::Two<T>(),nt2::Two<T>(),nt2::Three<T>()), nt2::Seven<T>(), 0);
  NT2_TEST_ULP_EQUAL(fma(nt2::Zero<T>(), nt2::Zero<T>(), nt2::Zero<T>()), nt2::Zero<T>(), 0);
} // end of test for signed_int_
