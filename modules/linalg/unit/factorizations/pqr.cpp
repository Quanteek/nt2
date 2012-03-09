//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II       
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI       
//                                                                              
//          Distributed under the Boost Software License, Version 1.0.          
//                 See accompanying file LICENSE.txt or copy at                 
//                     http://www.boost.org/LICENSE_1_0.txt                     
//==============================================================================
#define NT2_UNIT_MODULE "nt2 linalg toolbox - pqr"

//////////////////////////////////////////////////////////////////////////////
// unit test behavior of algebra components
//////////////////////////////////////////////////////////////////////////////

#include <nt2/table.hpp>
#include <nt2/include/functions/size.hpp>
#include <nt2/include/functions/pqr.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/make_vector.hpp>
#include <nt2/sdk/unit/tests.hpp>
#include <nt2/sdk/unit/module.hpp>

template < class TAB>
void pt(const TAB & t){}
NT2_TEST_CASE_TPL ( pqr, NT2_REAL_TYPES) 
{
  typedef nt2::table<T> table_t; 
  table_t b(nt2::of_size(4, 4));
  
  std::size_t dim1_b = size(b)(1);
  std::size_t dim2_b = size(b)(2);
  
  for(std::size_t i = 1; i <= dim1_b; i++)
    for(std::size_t j = 1; j <= dim2_b; j++)
      {
        b(i, j)= T((i == j)*i);
      }
  nt2::pqr_return<table_t> f = pqr(b);
//   std::cout << "1" << std::endl; 
//   table_t q = f.getq();
//   std::cout << "2" << std::endl; 
//   table_t r= f.getr();
//   std::cout << "3" << std::endl; 
//   for(std::size_t i = 1; i <= size(q, 1); i++)
//     {
//       for(std::size_t j = 1; j <= size(q, 2); j++)
//         {
//           std::cout << q(i, j) << ", "; 
//         }
//       std::cout << std::endl; 
//     }
//   std::cout << std::endl; 
//   for(std::size_t i = 1; i <= size(r, 1); i++)
//     {
//       for(std::size_t j = 1; j <= size(r, 2); j++)
//         {
//           std::cout << r(i, j) << ", "; 
//         }
//       std::cout << std::endl; 
//     }

}


