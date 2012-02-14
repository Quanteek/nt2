//////////////////////////////////////////////////////////////////////////////
///   Copyright 2003 and onward LASMEA UMR 6602 CNRS/U.B.P Clermont-Ferrand
///   Copyright 2009 and onward LRI    UMR 8623 CNRS/Univ Paris Sud XI
///
///          Distributed under the Boost Software License, Version 1.0
///                 See accompanying file LICENSE.txt or copy at
///                     http://www.boost.org/LICENSE_1_0.txt
//////////////////////////////////////////////////////////////////////////////
#define NT2_UNIT_MODULE "nt2 algebra toolbox - gemm"

//////////////////////////////////////////////////////////////////////////////
// timing Test behavior of blas gemm
//////////////////////////////////////////////////////////////////////////////
#include <nt2/toolbox/algebra/include/functions/gemm.hpp>
#include <nt2/sdk/memory/container.hpp>
#include <nt2/core/functions/of_size.hpp>
#include <nt2/include/functions/size.hpp>
#include <boost/simd/sdk/memory/align_on.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/make_vector.hpp>
#include <nt2/sdk/timing/now.hpp>
#include <nt2/sdk/unit/details/helpers.hpp>
#include <nt2/sdk/unit/perform_benchmark.hpp>
#include <nt2/sdk/unit/module.hpp>

using nt2::of_size;

template<class T> struct dispatch_gemm_test
{
  typedef nt2::memory::
          container <nt2::tag::table_, nt2::id_<0>, T, nt2::settings(nt2::_2D)>  table_2D_type;

  dispatch_gemm_test( long int m, long int n, T alpha, T beta) 
    : a(of_size(m,n)), b(of_size(n,m)), r(of_size(m,m))
    , AL(alpha), BE(beta) 
  {
    dim1_a = size(a)(1);
    dim2_a = size(a)(2);
    dim1_b = size(b)(1);
    dim2_b = size(b)(2);
    dim1_r = size(r)(1);
    dim2_r = size(r)(2);
    for(std::size_t i = 1; i <= dim1_a; i++)
      for(std::size_t j = 1; j <= dim2_a; j++)
        {
          pos1 = boost::fusion::make_vector(i,j);
          a[pos1] = i*j;
        }
    
    for(std::size_t i = 1; i <= dim1_b; i++)
      for(std::size_t j = 1; j <= dim2_b; j++)
        {
          pos1 = boost::fusion::make_vector(i,j);
          b[pos1] = i*j;
        }
    
    for(std::size_t i = 1; i <= dim1_a; i++)
      for(std::size_t j = 1; j <= dim2_b; j++)
        {
          pos1 = boost::fusion::make_vector(i,j);
          r[pos1] = 0;
        }
    
  }
  
  void operator()()
  {
    nt2::gemm(nt2::gemm_status<'N','N'>(), a, b, r, AL, BE);
  }

  table_2D_type a, b, r;
  boost::fusion::vector<int,int> pos1,pos2,pos3;
  T AL, BE;
  std::size_t dim1_a;
  std::size_t dim2_a;
  std::size_t dim1_b;
  std::size_t dim2_b;
  std::size_t dim1_r;
  std::size_t dim2_r;
  
};

template<class T> struct system_gemm_test
{
  typedef nt2::memory::
          container <nt2::tag::table_, nt2::id_<0>, T, nt2::settings(nt2::_2D)>  table_2D_type;
  typedef typename table_2D_type::parent::lead_t lead_t;

  system_gemm_test(  const char ta, const char tb, long int m
                   , long int n, T alpha, T beta) 
    : a(of_size(m,n)), b(of_size(n,m)), r(of_size(m,m))
    , M(nt2::size(a)(ta=='T'?2:1)), N(nt2::size(b)(tb=='T'?1:2)), K(nt2::size(a)(ta=='T'?1:2))
    , LDA(boost::simd::memory::align_on(size(a, 1), lead_t::value))
    , LDB(boost::simd::memory::align_on(size(b, 1), lead_t::value))
    , LDC(boost::simd::memory::align_on(size(a, 1), lead_t::value))
    , TA(ta), TB(tb), AL(alpha), BE(beta) 
  {
    dim1_a = size(a)(1);
    dim2_a = size(a)(2);
    dim1_b = size(b)(1);
    dim2_b = size(b)(2);
    dim1_r = size(r)(1);
    dim2_r = size(r)(2);
    for(std::size_t i = 1; i <= dim1_a; i++)
      for(std::size_t j = 1; j <= dim2_a; j++)
        {
          pos1 = boost::fusion::make_vector(i,j);
          a[pos1] = i*j;
        }
    
    for(std::size_t i = 1; i <= dim1_b; i++)
      for(std::size_t j = 1; j <= dim2_b; j++)
        {
          pos1 = boost::fusion::make_vector(i,j);
          b[pos1] = i*j;
        }
    
    for(std::size_t i = 1; i <= dim1_a; i++)
      for(std::size_t j = 1; j <= dim2_b; j++)
        {
          pos1 = boost::fusion::make_vector(i,j);
          r[pos1] = 0;
        }
    
  }
  
  void operator()()
  {
    nt2::ext::gemm(&TA, &TB, &M, &N, &K, &AL, a.begin(), &LDA, b.begin(), &LDB, &BE, r.begin(), &LDC);
  }

  table_2D_type a, b, r;
  long int M, N, K, LDA, LDB, LDC;
  boost::fusion::vector<int,int> pos1,pos2,pos3;
  T AL, BE;
  char TA, TB;
  std::size_t dim1_a;
  std::size_t dim2_a;
  std::size_t dim1_b;
  std::size_t dim2_b;
  std::size_t dim1_r;
  std::size_t dim2_r;
  
};

template<class T> struct normal_gemm_test
{
  typedef nt2::memory::
          container <nt2::tag::table_, nt2::id_<0>, T, nt2::settings(nt2::_2D)>  table_2D_type;
  typedef typename table_2D_type::parent::lead_t lead_t;

  normal_gemm_test(long int m, long int n) 
    : a(of_size(m,n)), b(of_size(n,m)), r(of_size(m,m))
  {
    dim1_a = size(a)(1);
    dim2_a = size(a)(2);
    dim1_b = size(b)(1);
    dim2_b = size(b)(2);
    dim1_r = size(r)(1);
    dim2_r = size(r)(2);
    for(std::size_t i = 1; i <= dim1_a; i++)
      for(std::size_t j = 1; j <= dim2_a; j++)
        {
          pos1 = boost::fusion::make_vector(i,j);
          a[pos1] = i*j;
        }
    
    for(std::size_t i = 1; i <= dim1_b; i++)
      for(std::size_t j = 1; j <= dim2_b; j++)
        {
          pos1 = boost::fusion::make_vector(i,j);
          b[pos1] = i*j;
        }
    
    for(std::size_t i = 1; i <= dim1_a; i++)
      for(std::size_t j = 1; j <= dim2_b; j++)
        {
          pos1 = boost::fusion::make_vector(i,j);
          r[pos1] = 0;
        }
    
  }
  
  void operator()()
  {
    T tmp;
    for(std::size_t i = 1; i <= dim1_a; i++)
      for(std::size_t k = 1; k <= dim2_a; k++)
      {
        pos1 = boost::fusion::make_vector(i,k);
        tmp  = a[pos1];
        for(std::size_t j = 1; j <= dim2_b; j++)
        {
          pos2 = boost::fusion::make_vector(i,j);
          pos3 = boost::fusion::make_vector(k,j);
          r[pos2] += tmp*b[pos3];
        }
      }
  }

  table_2D_type a, b, r;
  boost::fusion::vector<int,int> pos1,pos2,pos3;
  std::size_t dim1_a;
  std::size_t dim2_a;
  std::size_t dim1_b;
  std::size_t dim2_b;
  std::size_t dim1_r;
  std::size_t dim2_r;
  
};

template<class T> void do_test()
{
  std::cout << "============================\n";
  std::cout << "Matrix product dummy\n";
  std::cout << "============================\n";
  normal_gemm_test<T> ngt(1024, 512);
  double dv1 = nt2::unit::perform_benchmark( ngt, 1.);
  std::cout << dv1/(1024*1024) << "\t\n";
 
  std::cout << "============================\n";
  std::cout << "Dispatch blas\n";
  std::cout << "============================\n";
  dispatch_gemm_test<T> dgt(1024, 512, 1.0, 0.0);
  double dv3 = nt2::unit::perform_benchmark( dgt, 1.);
  std::cout << dv3/(1024*1024) << "\t\n";

  std::cout << "============================\n";
  std::cout << "System blas\n";
  std::cout << "============================\n";
  system_gemm_test<T> sgt('N','N', 1024, 512, 1.0, 0.0);
  double dv2 = nt2::unit::perform_benchmark( sgt, 1.);
  std::cout << dv2/(1024*1024) << "\t\n";
}

NT2_TEST_CASE_TPL( gemm, (double)(float) )
{
  do_test<T>();
}
