//////////////////////////////////////////////////////////////////////////////
///   Copyright 2003 and onward LASMEA UMR 6602 CNRS/U.B.P Clermont-Ferrand
///   Copyright 2009 and onward LRI    UMR 8623 CNRS/Univ Paris Sud XI
///
///          Distributed under the Boost Software License, Version 1.0
///                 See accompanying file LICENSE.txt or copy at
///                     http://www.boost.org/LICENSE_1_0.txt
//////////////////////////////////////////////////////////////////////////////
#define NT2_UNIT_MODULE "nt2 trigonometric toolbox - unit/simd Mode"

#include <nt2/toolbox/trigonometric/include/fast_cospi.hpp>
#include <nt2/sdk/unit/tests.hpp>
#include <nt2/sdk/unit/module.hpp>
#include <nt2/sdk/simd/native.hpp>
#include <nt2/sdk/memory/is_aligned.hpp>
#include <nt2/sdk/memory/aligned_type.hpp>
#include <nt2/sdk/memory/load.hpp>
#include <nt2/sdk/functor/meta/call.hpp>
#include <boost/type_traits/is_same.hpp>
#include <nt2/include/functions/random.hpp>
#include <nt2/include/functions/ulpdist.hpp>
#include <iostream>

//////////////////////////////////////////////////////////////////////////////
// Test behavior of arithmetic components using NT2_TEST_CASE
//////////////////////////////////////////////////////////////////////////////
NT2_TEST_CASE_TPL(fast_cospi, (double)(float))//NT2_SIMD_REAL_CONVERTIBLE_TYPES )
{
 using nt2::fast_cospi;
 using nt2::functors::fast_cospi_;    
 using nt2::load; 
 using nt2::simd::native; 
 using nt2::meta::cardinal_of;

 typedef NT2_SIMD_DEFAULT_EXTENSION  ext_t;
 typedef native<T,ext_t>             n_t;
 typedef typename nt2::meta::call<fast_cospi_(n_t)>::type call_type;
 typedef typename nt2::meta::as_real<T>::type rT;
 typedef native<rT,ext_t>             rn_t;
 
   
 NT2_TEST( (boost::is_same<call_type, rn_t>::value) );  
 NT2_ALIGNED_TYPE(T) data[1*cardinal_of<n_t>::value];
 double z, m = 0; 
 for(int num = 0; num < 10; num++)
   {
     for(std::size_t i=0;i<1*cardinal_of<n_t>::value;++i){
       data[i] = nt2::random(0.0, 0.25); // good value here for fast_cospi
     }
     n_t a0 = load<n_t>(&data[0],0); 
     rn_t v  = fast_cospi(a0);
     for(std::size_t j=0;j<cardinal_of<n_t>::value;++j)
       {
	 NT2_TEST_LESSER( z = nt2::ulpdist(v[j], fast_cospi(a0[j])), 1);
	 if (z > m) m = z; 
       }
   }
 std::cout << "ulp max = " << m << std::endl;
}
