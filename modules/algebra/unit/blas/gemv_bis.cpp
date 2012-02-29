//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II       
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI       
//                                                                              
//          Distributed under the Boost Software License, Version 1.0.          
//                 See accompanying file LICENSE.txt or copy at                 
//                     http://www.boost.org/LICENSE_1_0.txt                     
//==============================================================================
#define NT2_UNIT_MODULE "nt2 algebra toolbox - gemv"

//////////////////////////////////////////////////////////////////////////////
// unit test behavior of algebra components
//////////////////////////////////////////////////////////////////////////////

#include <nt2/table.hpp>
#include <nt2/include/functions/gemv.hpp>
#include <nt2/include/functions/size.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/make_vector.hpp>
#include <nt2/sdk/unit/tests.hpp>
#include <nt2/sdk/unit/module.hpp>
#include <nt2/toolbox/algebra/details/padding.hpp>

NT2_TEST_CASE_TPL ( gemv, NT2_REAL_TYPES) 
{
  using nt2::meta::make_container; 
  using nt2::of_size_;
  using nt2::gemv;
  typedef typename make_container<nt2::tag::table_, T, of_size_<7, 5>  >::type table_type_a;
  typedef typename make_container<nt2::tag::table_, T, of_size_<1, 5> >::type table_type_b;
  typedef typename make_container<nt2::tag::table_, T, of_size_<7> >::type table_type_r;

  table_type_a a;
  table_type_b b;
  table_type_r r;
  table_type_r r_;
  boost::fusion::vector<int,int> pos1;
  boost::fusion::vector<int,int> pos2;
  boost::fusion::vector<int> pos3;
  std::size_t dim1_a = size(a)(1);
  std::size_t dim2_a = size(a)(2);
  T tmp;

  for(std::size_t i = 1; i <= dim1_a; i++)
    for(std::size_t j = 1; j <= dim2_a; j++)
    {
      pos1 = boost::fusion::make_vector(i,j);
      a[pos1] =i*j;
    }

  for(std::size_t j = 1; j <= dim2_a; j++)
  {
    pos2 = boost::fusion::make_vector(1, j);
    b[pos2] = j;
  }

  for(std::size_t j = 1; j <= dim2_a; j++)
  {
    pos3 = boost::fusion::make_vector(j);
    r[pos3] = r_[pos3] = 0.0;
  }

  for(std::size_t i = 1; i <= dim1_a; i++)
    {
      pos3 = boost::fusion::make_vector(i);
      for(std::size_t k = 1; k <= dim2_a; k++)
        {
          pos2 = boost::fusion::make_vector(1, k);
          pos1 = boost::fusion::make_vector(i, k);
          r[pos3] += a[pos1]*b[pos2];
        }
    }
  
  // Call blas
  gemv(nt2::gemv_status<'N', 'T'>(), a, b, r_);

  for(std::size_t j = 1; j <= dim1_a; j++)
  {
    pos3 = boost::fusion::make_vector(j);
    NT2_TEST_EQUAL(r[pos3], r_[pos3]);
  }
//   const T * z =  b.begin();
//   for(int i = 0; i < 5; ++i)
//     {
//     std::cout << z << "  " << *z << "  "  << std::endl;
//     z+= nt2::details::padding(b); 
//     }
//   std::cout << std::endl; 
//   for(int i = 1; i <= 5; ++i){
//     pos2 = boost::fusion::make_vector(1, i);
//     std::cout << "  " << b[pos2];
//     std::cout << "  " << &b[pos2];
//   }
//   std::cout << std::endl<< std::endl; 
//   std::cout << std::endl; 
//   for(int i = 2; i <= 5; ++i){
//     pos3 = boost::fusion::make_vector(1, i);
//     pos2 = boost::fusion::make_vector(1, i-1);
//     std::cout << "  " << &b[pos2]<< "  " << &b[pos3] << std::endl; 
//     std::cout << "  " << ((long int)(&b[pos3]) - (long int)(&b[pos2]))<< std::endl; 
//   }
}
