//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II       
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI       
//                                                                              
//          Distributed under the Boost Software License, Version 1.0.          
//                 See accompanying file LICENSE.txt or copy at                 
//                     http://www.boost.org/LICENSE_1_0.txt                     
//==============================================================================
#define NT2_UNIT_MODULE "nt2 linalg toolbox - chol_ip"

//////////////////////////////////////////////////////////////////////////////
// unit test behavior of algebra components
//////////////////////////////////////////////////////////////////////////////

#include <nt2/table.hpp>
#include <nt2/include/functions/size.hpp>
#include <nt2/toolbox/linalg/functions/solvers/solve_chol_ip.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/make_vector.hpp>
#include <nt2/sdk/unit/tests.hpp>
#include <nt2/sdk/unit/module.hpp>

template < class TAB>
void pt(const TAB & t){}
NT2_TEST_CASE_TPL ( LUSolveIP, (double)) 
{
  typedef nt2::table<T> table_t; 
  table_t a(nt2::of_size(4, 4));
  table_t x(nt2::of_size(4, 1));
  
  for(std::size_t i = 1; i <= size(a, 1); i++)
    {
      x(i, 1) = T(i); 
      for(std::size_t j = 1; j <= size(a, 2); j++)
        {
          a(i, j)= T((i == j)*i);
        }
    }
  a(1, 4) = T(1);
  a(4, 1) = T(1);
  x(1) = T(2);
  x(4) = T(5); 
  for(std::size_t i = 1; i <= size(a, 1); i++)
    {
      for(std::size_t j = 1; j <= size(a, 2); j++)
        {
          std::cout << a(i, j)<< ", "; 
        }
      std::cout << std::endl;
    }
  std::cout << std::endl;
  for(std::size_t i = 1; i <= size(x, 1); i++)
    {
      std::cout << x(i, 1) << ", "; 
    }
  std::cout << std::endl;
  nt2::solve_chol_ip_return<table_t, table_t> f = solve_chol_ip(a, x);
  std::cout <<  "a ->"<< std::endl; 
  for(std::size_t i = 1; i <= size(a, 1); i++)
    {
      for(std::size_t j = 1; j <= size(a, 2); j++)
        {
          std::cout << a(i, j)<< ", "; 
        }
      std::cout << std::endl;
    }
  std::cout << std::endl;
  std::cout <<  "x -> "; 
  for(std::size_t i = 1; i <= size(x, 1); i++)
    {
      std::cout << x(i) << ", "; 
    }
}


