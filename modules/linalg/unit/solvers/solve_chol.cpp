//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II       
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI       
//                                                                              
//          Distributed under the Boost Software License, Version 1.0.          
//                 See accompanying file LICENSE.txt or copy at                 
//                     http://www.boost.org/LICENSE_1_0.txt                     
//==============================================================================
#define NT2_UNIT_MODULE "nt2 linalg toolbox - solve_chol"

//////////////////////////////////////////////////////////////////////////////
// unit test behavior of algebra components
//////////////////////////////////////////////////////////////////////////////

#include <nt2/table.hpp>
#include <nt2/include/functions/size.hpp>
#include <nt2/toolbox/linalg/functions/solvers/solve.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/make_vector.hpp>
#include <nt2/sdk/unit/tests.hpp>
#include <nt2/sdk/unit/module.hpp>

NT2_TEST_CASE_TPL ( CHOLSolve, (double)) 
{
  typedef nt2::table<T> table_t; 
  table_t a(nt2::of_size(4, 4));
  table_t b(nt2::of_size(4, 1));
  
  for(std::size_t i = 1; i <= size(a, 1); i++)
    {
      b(i, 1) = T(i); 
      for(std::size_t j = 1; j <= size(a, 2); j++)
        {
          a(i, j)= T((i == j)*i);
        }
    }
  a(1, 4) = T(1);
  a(4, 1) = T(1);
  b(1) = T(2);
  b(4) = T(5); 
  for(std::size_t i = 1; i <= size(a, 1); i++)
    {
      for(std::size_t j = 1; j <= size(a, 2); j++)
        {
          std::cout << a(i, j)<< ", "; 
        }
      std::cout << std::endl;
    }
  std::cout << std::endl;
  for(std::size_t i = 1; i <= size(b, 1); i++)
    {
      std::cout << b(i) << ", "; 
    }
  std::cout << std::endl;
  std::cout << "saving data" << std::endl;
  //saving data
  table_t x = solve_chol(a, b);
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
  std::cout <<  "b -> "; 
  for(std::size_t i = 1; i <= size(b, 1); i++)
    {
      std::cout << b(i) << ", "; 
    }
  std::cout <<  "x -> "; 
  for(std::size_t i = 1; i <= size(x, 1); i++)
    {
      std::cout << x(i) << ", "; 
    }
  std::cout << std::endl;
  std::cout << "destroying data" << std::endl;
  //destroying data  
  x = solve_chol(a, b, nt2::allowdestroy());
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
  std::cout <<  "b -> "; 
  for(std::size_t i = 1; i <= size(b, 1); i++)
    {
      std::cout << b(i) << ", "; 
    }
  std::cout <<  "x -> "; 
  for(std::size_t i = 1; i <= size(x, 1); i++)
    {
      std::cout << x(i) << ", "; 
    }  
}


