//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II       
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI       
//                                                                              
//          Distributed under the Boost Software License, Version 1.0.          
//                 See accompanying file LICENSE.txt or copy at                 
//                     http://www.boost.org/LICENSE_1_0.txt                     
//==============================================================================
#define NT2_UNIT_MODULE "nt2 linalg toolbox - chol"

//////////////////////////////////////////////////////////////////////////////
// unit test behavior of algebra components
//////////////////////////////////////////////////////////////////////////////

#include <nt2/table.hpp>
#include <nt2/include/functions/size.hpp>
#include <nt2/include/functions/chol.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/make_vector.hpp>
#include <nt2/sdk/unit/tests.hpp>
#include <nt2/sdk/unit/module.hpp>
#include <iostream>

template < class TAB>
void pt(const TAB & t){}
NT2_TEST_CASE_TPL ( chol, NT2_REAL_TYPES) 
{
  typedef nt2::table<T> table_t; 
  table_t b(nt2::of_size(4, 4));
  
  std::size_t dim1_b = size(b)(1);
  std::size_t dim2_b = size(b)(2);
  
  for(std::size_t i = 1; i <= dim1_b; i++)
    for(std::size_t j = 1; j <= dim2_b; j++)
      {
        b(i, j)= T(1+10*(i == j));
      }
  nt2::chol_f<table_t> f = chol(b);
  table_t u = f.getu();
  for(std::size_t i = 1; i <= size(u, 1); i++)
    {
      for(std::size_t j = 1; j <= size(u, 2); j++)
        {
          std::cout << u(i, j) << "\t";
        }
      std::cout  << std::endl; 
    }
  std::cout << "rcond =  "   << f.rcond() << std::endl;
  std::cout << "negativity " << f.negativity() << std::endl;
}

                                 
