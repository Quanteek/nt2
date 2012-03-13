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
#include <nt2/toolbox/linalg/functions/solvers/internal.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/make_vector.hpp>
#include <nt2/sdk/unit/tests.hpp>
#include <nt2/sdk/unit/module.hpp>

template < class TAB>
void pt(const TAB & t){}
NT2_TEST_CASE_TPL ( LUSolveIP, NT2_REAL_TYPES) 
{
  typedef nt2::table<T> table_t; 
  table_t b(nt2::of_size(4, 4));
  table_t x(nt2::of_size(4, 1));
  table_t xx(nt2::of_size(4, 1));
  
  for(std::size_t i = 1; i <= size(b, 1); i++)
    {
      x(i, 1) = T(i); 
      for(std::size_t j = 1; j <= size(b, 2); j++)
        {
          b(i, j)= T((i == j)*i);
        }
    }
  for(std::size_t i = 1; i <= size(b, 1); i++)
    {
      for(std::size_t j = 1; j <= size(b, 2); j++)
        {
          std::cout << b(i, j)<< ", "; 
        }
      std::cout << std::endl;
    }
  std::cout << std::endl;
  for(std::size_t i = 1; i <= size(x, 1); i++)
    {
      std::cout << x(i, 1) << ", "; 
    }
  std::cout << std::endl;
  nt2::details::internal_solvers< table_t, table_t>::QRSolveIP(b, xx, x); 
  for(std::size_t i = 1; i <= size(b, 1); i++)
    {
      for(std::size_t j = 1; j <= size(b, 2); j++)
        {
          std::cout << b(i, j)<< ", "; 
        }
      std::cout << std::endl;
    }
  std::cout << std::endl;
  for(std::size_t i = 1; i <= size(xx, 1); i++)
    {
      std::cout << xx(i, 1) << ", "; 
    }
  std::cout << std::endl; 
}


