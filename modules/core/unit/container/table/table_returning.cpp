/*******************************************************************************
 *         Copyright 2003 & onward LASMEA UMR 6602 CNRS/Univ. Clermont II
 *         Copyright 2009 & onward LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#define NT2_UNIT_MODULE "nt2::table_returning "

#include <nt2/table.hpp>
#include <nt2/include/functions/ulpdist.hpp>
#include <nt2/include/functions/first_index.hpp>
#include <boost/type_traits/common_type.hpp>
#include <nt2/sdk/unit/tests.hpp>
#include <nt2/sdk/unit/module.hpp>

////////////////////////////////////////////////////////////////////////////////
// table type has some dimensions
////////////////////////////////////////////////////////////////////////////////

template < class T > T transmit(const T& a)
{
  return a; 
}
NT2_TEST_CASE_TPL( table_dimensions ,NT2_TYPES)
{
  using nt2::table;
  double ulpd = 0.0;
  typedef T r_t; 
  nt2::table<T> a(nt2::of_size(3, 3));
  
 for(int i=1; i <= 3; i++)
   {
     for(int j=1; j <= 3; j++)
       {
         a(i, j) = T(i + j); 
       }
   }
 nt2::table<T, nt2::C_index_> b(nt2::of_size(3, 3)); 
 //  b = transmit(a);
  b =  a;
  int fa1 = nt2::first_index<1>(a);
  int fb1 = nt2::first_index<1>(b); 
  int fa2 = nt2::first_index<2>(a);
  int fb2 = nt2::first_index<2>(b); 
  for(int i=0; i < 3; i++)
    {
      for(int j=0; j < 3; j++)
        {
          NT2_TEST_ULP_EQUAL(a(i+fa1, j+fa2), b(i+fb1, j+fb2), 0); 
        }
    }
  
  
}



