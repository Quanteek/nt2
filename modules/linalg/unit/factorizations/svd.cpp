//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II       
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI       
//                                                                              
//          Distributed under the Boost Software License, Version 1.0.          
//                 See accompanying file LICENSE.txt or copy at                 
//                     http://www.boost.org/LICENSE_1_0.txt                     
//==============================================================================
#define NT2_UNIT_MODULE "nt2 linalg toolbox - svd"

//////////////////////////////////////////////////////////////////////////////
// unit test behavior of algebra components
//////////////////////////////////////////////////////////////////////////////

#include <nt2/table.hpp>
#include <nt2/include/functions/size.hpp>
#include <nt2/include/functions/svd.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/make_vector.hpp>
#include <nt2/sdk/unit/tests.hpp>
#include <nt2/sdk/unit/module.hpp>

template < class TAB>
void pt(const TAB & t){}
NT2_TEST_CASE_TPL ( svd, NT2_REAL_TYPES) 
{
  typedef nt2::table<T> table_t; 
  table_t b(nt2::of_size(4, 4));
  
  std::size_t dim1_b = size(b, 1);
  std::size_t dim2_b = size(b, 2);
  
  for(std::size_t i = 1; i <= dim1_b; i++)
    for(std::size_t j = 1; j <= dim2_b; j++)
      {
        b(i, j)= T((i == j)*i);
      }
  b(1, dim2_b) = T(1); 
  for(std::size_t i = 1; i <= size(b, 1); i++)
    {
      for(std::size_t j = 1; j <= size(b, 2); j++)
        {
          std::cout << b(i, j) << ", "; 
        }
      std::cout << std::endl; 
    }
  nt2::svd_return<table_t> f = svd(b);
  table_t u = f.get_u();
  table_t vt= f.get_vt();
  table_t w = f.get_singular(); 
  for(std::size_t i = 1; i <= size(u, 1); i++)
    {
      for(std::size_t j = 1; j <= size(u, 2); j++)
        {
          std::cout << u(i, j) << ", "; 
        }
      std::cout << std::endl; 
    }
  std::cout << std::endl; 
  for(std::size_t i = 1; i <= size(vt, 1); i++)
    {
      for(std::size_t j = 1; j <= size(vt, 2); j++)
        {
          std::cout << vt(i, j) << ", "; 
        }
      std::cout << std::endl; 
    }
  std::cout << std::endl; 
  for(std::size_t i = 1; i <= nt2::numel(w); i++)
    std::cout << w(i) << ", "; 
}


