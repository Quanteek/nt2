//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II       
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI       
//                                                                              
//          Distributed under the Boost Software License, Version 1.0.          
//                 See accompanying file LICENSE.txt or copy at                 
//                     http://www.boost.org/LICENSE_1_0.txt                     
//==============================================================================
#define NT2_UNIT_MODULE "nt2 linalg toolbox - lse"

//////////////////////////////////////////////////////////////////////////////
// unit test behavior of algebra components
//////////////////////////////////////////////////////////////////////////////

#include <nt2/table.hpp>
#include <nt2/include/functions/lse.hpp>
#include <nt2/include/functions/size.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/make_vector.hpp>
#include <nt2/sdk/unit/tests.hpp>
#include <nt2/sdk/unit/module.hpp>

NT2_TEST_CASE_TPL ( lse, NT2_REAL_TYPES) 
{
  using nt2::meta::make_container; 
  using nt2::of_size_;
  using nt2::lse;
  const int N = 5;
  const int M = 5;
  const int P = 7; 
  typedef typename make_container<nt2::tag::table_, T, of_size_<M, M>  >::type table_type_a;
  typedef typename make_container<nt2::tag::table_, T, of_size_<M, 2>  >::type table_type_b;
  typedef typename make_container<nt2::tag::table_, T, of_size_<M, 2>  >::type table_type_c;
  typedef typename make_container<nt2::tag::table_, T, of_size_<M, 2>  >::type table_type_d;
  typedef typename make_container<nt2::tag::table_, T, of_size_<M, 2>  >::type table_type_x;
  table_type_a a;
  table_type_b b;
  table_type_c c;
  table_type_x x;
  table_type_d d; 
  boost::fusion::vector<int,int> pos1,pos2,pos3;
  std::size_t dim1_a = size(a)(1);
  std::size_t dim2_a = size(a)(2);
  std::size_t dim1_b = size(b)(1);
  std::size_t dim2_b = size(b)(2);
//   std::size_t dim1_r = size(r)(1);
//   std::size_t dim2_r = size(r)(2);
//   std::size_t dim1_r_ = size(r)(1);
//   std::size_t dim2_r_ = size(r)(2);
//   T tmp;

//   for(std::size_t i = 1; i <= dim1_a; i++)
//     for(std::size_t j = 1; j <= dim2_a; j++)
//     {
//       pos1 = boost::fusion::make_vector(i,j);
//       a[pos1] = i*j;
//     }

//   for(std::size_t i = 1; i <= dim1_b; i++)
//     for(std::size_t j = 1; j <= dim2_b; j++)
//     {
//       pos1 = boost::fusion::make_vector(i,j);
//       b[pos1] = i*j;
//     }

//   for(std::size_t i = 1; i <= dim1_a; i++)
//     for(std::size_t j = 1; j <= dim2_b; j++)
//     {
//       pos1 = boost::fusion::make_vector(i,j);
//       r[pos1] = 0;
//     }

  

//   for(std::size_t i = 1; i <= dim1_a; i++)
//     for(std::size_t k = 1; k <= dim2_a; k++)
//     {
//       pos1 = boost::fusion::make_vector(i,k);
//       tmp  = a[pos1];
//       for(std::size_t j = 1; j <= dim2_b; j++)
//       {
//         pos2 = boost::fusion::make_vector(i,j);
//         pos3 = boost::fusion::make_vector(k,j);
//         r[pos2] += tmp*b[pos3];
//       }
//     }

//   {
//     for(std::size_t i = 1; i <= dim1_r_; i++)
//       for(std::size_t j = 1; j <= dim2_r_; j++)
//         {
//           pos1 = boost::fusion::make_vector(i,j);
//           r_[pos1] = 0;
//         }
    
//     // Call blas
//     nt2::lse(nt2::b_status<>(), a, b, r_);
    
//     for(std::size_t i = 1; i <= dim1_r; i++)
//       for(std::size_t j = 1; j <= dim2_r; j++)
//         {
//           pos1 = boost::fusion::make_vector(i,j);
//           std::cout << "i " << i << " j " << j << " -> " << r[pos1] << ", " << r_[pos1] << std::endl; 
//           NT2_TEST_EQUAL(r[pos1], r_[pos1]);
//         }
//   }
//   {
//     for(std::size_t i = 1; i <= dim1_r_; i++)
//       for(std::size_t j = 1; j <= dim2_r_; j++)
//         {
//           pos1 = boost::fusion::make_vector(i,j);
//           r_[pos1] = 0;
//         }
    std::cout << " ==================================================== " << std::endl; 
    // Call blas
    nt2::lse(a, b, c, d, x);
    
//     for(std::size_t i = 1; i <= dim1_r; i++)
//       for(std::size_t j = 1; j <= dim2_r; j++)
//         {
//           pos1 = boost::fusion::make_vector(i,j);
//           std::cout << "i " << i << " j " << j << " -> " << r[pos1] << ", " << r_[pos1] << std::endl; 
//           NT2_TEST_EQUAL(r[pos1], r_[pos1]);
//         }
//   }

  
}

