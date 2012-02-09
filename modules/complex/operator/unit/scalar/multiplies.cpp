//////////////////////////////////////////////////////////////////////////////
///   Copyright 2003 and onward LASMEA UMR 6602 CNRS/U.B.P Clermont-Ferrand
///   Copyright 2009 and onward LRI    UMR 8623 CNRS/Univ Paris Sud XI
///
///          Distributed under the Boost Software License, Version 1.0
///                 See accompanying file LICENSE.txt or copy at
///                     http://www.boost.org/LICENSE_1_0.txt
//////////////////////////////////////////////////////////////////////////////
#define NT2_UNIT_MODULE "nt2 complex.operator toolbox - multiplies/scalar Mode"

//////////////////////////////////////////////////////////////////////////////
// unit test behavior of boost.simd.operator components in scalar mode
//////////////////////////////////////////////////////////////////////////////
/// created  by jt the 18/02/2011
/// 
#include <nt2/include/functions/multiplies.hpp>
#include <nt2/include/functions/ulpdist.hpp>
#include <boost/simd/sdk/simd/logical.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/dispatch/functor/meta/call.hpp>
#include <nt2/sdk/unit/tests.hpp>
#include <nt2/sdk/unit/module.hpp>
#include <boost/simd/sdk/memory/buffer.hpp>
#include <nt2/toolbox/constant/constant.hpp>

NT2_TEST_CASE_TPL ( multiplies_real__2_0,  BOOST_SIMD_REAL_TYPES)
{
  
  using nt2::multiplies;
  using nt2::tag::multiplies_;
  typedef std::complex<T> cT; 
  typedef typename boost::dispatch::meta::as_integer<T>::type iT;
  typedef typename boost::dispatch::meta::call<multiplies_(cT,cT)>::type r_t;
  typedef typename nt2::meta::scalar_of<r_t>::type sr_t;
  typedef typename nt2::meta::scalar_of<r_t>::type ssr_t;
  typedef typename boost::dispatch::meta::upgrade<T>::type u_t;
  typedef nt2::imaginary<T> ciT; 
  typedef cT wished_r_t;

  // return type conformity test 
  NT2_TEST( (boost::is_same < r_t, wished_r_t >::value) );
  std::cout << std::endl; 
  double ulpd;
  ulpd=0.0;


  // specific values tests
//   cT infz = cT(nt2::Inf<T>());
//   cT zinf = cT(0, nt2::Inf<T>());  
//   cT nanz = cT(nt2::Nan<T>());
//   cT znan = cT(0, nt2::Nan<T>());
//   cT z    = cT(0, 0); 
//   NT2_TEST_EQUAL(nt2::multiplies(infz,zinf), zinf);
//   std::cout << nt2::multiplies(infz,zinf) << std::endl; 
//   NT2_TEST_EQUAL(nt2::multiplies(infz,z), nanz);
//   std::cout << nt2::multiplies(infz, z) << std::endl; 
//   NT2_TEST_EQUAL(nt2::multiplies(zinf, z), znan);
//   std::cout << nt2::multiplies(zinf, z) << std::endl; 
//   NT2_TEST_EQUAL(nt2::multiplies(cT(nt2::Inf<T>()),  cT(nt2::Inf<T>())),  cT(nt2::Inf<T>()));
//   cT infinf = cT(nt2::Inf<T>(),nt2::Inf<T>()); 
//   NT2_TEST_EQUAL(nt2::multiplies(infinf,infinf),  cT(nt2::Nan<T>(), nt2::Inf<T>()));
//   std::cout << nt2::multiplies(infinf,infinf) << std::endl; 
//   NT2_TEST_EQUAL(nt2::multiplies(cT(nt2::Inf<T>()),  cT(nt2::Inf<T>())),  cT(nt2::Inf<T>()));
//   std::cout << nt2::multiplies(cT(nt2::Inf<T>()),  cT(nt2::Inf<T>())) << std::endl; 
//   NT2_TEST_EQUAL(nt2::multiplies(cT(nt2::Minf<T>()), cT(nt2::Minf<T>())), cT(nt2::Inf<T>()));
//   std::cout << nt2::multiplies(cT(nt2::Minf<T>()), cT(nt2::Minf<T>()))<< std::endl; 
//   NT2_TEST_EQUAL(nt2::multiplies(cT(nt2::Nan<T>()),  cT(nt2::Nan<T>())),  cT(nt2::Nan <T>()));   
//   std::cout <<nt2::multiplies(cT(nt2::Nan<T>()),  cT(nt2::Nan<T>())) << std::endl; 
//   NT2_TEST_EQUAL(nt2::multiplies(cT(nt2::One<T>()),  cT(nt2::Zero<T>())), cT(nt2::Zero<T>())); 
//   std::cout <<nt2::multiplies(cT(nt2::One<T>()),  cT(nt2::Zero<T>())) << std::endl; 
//   NT2_TEST_EQUAL(nt2::multiplies(cT(nt2::Zero<T>()), cT(nt2::Zero<T>())), cT(nt2::Zero<T>())); 
//   std::cout << nt2::multiplies(cT(nt2::Zero<T>()), cT(nt2::Zero<T>())) << std::endl; 
//   NT2_TEST_EQUAL(nt2::multiplies(cT(0, 1), cT(0, 1)), cT(-1, 0));
//   NT2_TEST_EQUAL(nt2::multiplies(cT(1, 0), T(1))    , cT(1, 0)); 
//   NT2_TEST_EQUAL(nt2::multiplies(cT(0, 2), cT(0, 1)), cT(-2, 0));
//   NT2_TEST_EQUAL(nt2::multiplies(cT(0, 1), ciT(1))   ,cT(-1, 0)); 
//   NT2_TEST_EQUAL(nt2::multiplies(ciT(1), ciT(1))     , T(-1)); 
//   NT2_TEST_EQUAL(nt2::multiplies(ciT(0), ciT(1))     , T(0));
//   NT2_TEST_EQUAL(nt2::multiplies(T(1),   ciT(2))     , ciT(2));

  NT2_TEST_EQUAL(nt2::multiplies(cT(nt2::Inf<T>(),  nt2::Zero<T>()), cT(nt2::Zero<T>(), nt2::Zero<T>())), cT(nt2::Nan<T>(),  nt2::Zero<T>())); 
  NT2_TEST_EQUAL(nt2::multiplies(cT(nt2::Inf<T>(),  nt2::Zero<T>()), cT(nt2::Zero<T>(), nt2::Inf<T>())),  cT(nt2::Zero<T>(), nt2::Inf<T>() )); 
  NT2_TEST_EQUAL(nt2::multiplies(cT(nt2::Inf<T>(),  nt2::Zero<T>()), cT(nt2::One<T>(),  nt2::Zero<T>())), cT(nt2::Inf<T>(),  nt2::Zero<T>())); 
  NT2_TEST_EQUAL(nt2::multiplies(cT(nt2::Inf<T>(),  nt2::Zero<T>()), cT(nt2::Inf<T>(),  nt2::Zero<T>())), cT(nt2::Inf<T>(),  nt2::Zero<T>())); 
  NT2_TEST_EQUAL(nt2::multiplies(cT(nt2::Inf<T>(),  nt2::Inf<T>()),  cT(nt2::Inf<T>(),  nt2::Zero<T>())), cT(nt2::Inf<T>(),  nt2::Inf<T>() )); 
  NT2_TEST_EQUAL(nt2::multiplies(cT(nt2::Inf<T>(),  nt2::Zero<T>()), cT(nt2::Inf<T>(),  nt2::Inf<T>())),  cT(nt2::Inf<T>(),  nt2::Inf<T>() )); 
  NT2_TEST_EQUAL(nt2::multiplies(cT(nt2::Inf<T>(),  nt2::Inf<T>()),  cT(nt2::Inf<T>(),  nt2::Inf<T>())),  cT(nt2::Nan<T>(),  nt2::Inf<T>() )); 
  NT2_TEST_EQUAL(nt2::multiplies(cT(nt2::Inf<T>(),  nt2::Zero<T>()), cT(nt2::One<T>(),  nt2::Zero<T>())), cT(nt2::Inf<T>(),  nt2::Zero<T>())); 
  NT2_TEST_EQUAL(nt2::multiplies(cT(nt2::Inf<T>(),  nt2::Zero<T>()), cT(nt2::One<T>(),  nt2::Inf<T>())),  cT(nt2::Inf<T>(),  nt2::Inf<T>() )); 
  NT2_TEST_EQUAL(nt2::multiplies(cT(nt2::Inf<T>(),  nt2::Zero<T>()), cT(nt2::One<T>(),  nt2::One<T>())),  cT(nt2::Inf<T>(),  nt2::Inf<T>() )); 
  NT2_TEST_EQUAL(nt2::multiplies(cT(nt2::Inf<T>(),  nt2::Inf<T>()),  cT(nt2::One<T>(),  nt2::Zero<T>())), cT(nt2::Inf<T>(),  nt2::Inf<T>() ));
  NT2_TEST_EQUAL(nt2::multiplies(cT(nt2::One<T>(),  nt2::Zero<T>()), cT(nt2::Inf<T>(),  nt2::Inf<T>())) , cT(nt2::Inf<T>(),  nt2::Inf<T>() ));
  NT2_TEST_EQUAL(nt2::multiplies(cT(nt2::Inf<T>(),  nt2::Inf<T>()),  cT(nt2::One<T>(),  nt2::One<T>())),  cT(nt2::Nan<T>(),  nt2::Inf<T>() )); 
  NT2_TEST_EQUAL(nt2::multiplies(cT(nt2::One<T>(),  nt2::Inf<T>()),  cT(nt2::Zero<T>(),  nt2::One<T>())),  cT(nt2::Minf<T>(),  nt2::One<T>() )); 
  NT2_TEST_EQUAL(nt2::multiplies(cT(nt2::Inf<T>(),  nt2::One<T>()),  cT(nt2::Zero<T>(),  nt2::One<T>())),  cT(nt2::Mone<T>(),  nt2::Inf<T>() )); 
  std::cout << " std::complex<float>(nt2::One<T>(),  nt2::Inf<T>()) " << std::complex<float>(nt2::One<T>(),  nt2::Inf<T>()) << std::endl; 
 
} // end of test for floating_
