//////////////////////////////////////////////////////////////////////////////
///   Copyright 2003 and onward LASMEA UMR 6602 CNRS/U.B.P Clermont-Ferrand
///   Copyright 2009 and onward LRI    UMR 8623 CNRS/Univ Paris Sud XI
///
///          Distributed under the Boost Software License, Version 1.0
///                 See accompanying file LICENSE.txt or copy at
///                     http://www.boost.org/LICENSE_1_0.txt
//////////////////////////////////////////////////////////////////////////////
#define NT2_UNIT_MODULE "nt2 algebra toolbox - gemv"

//////////////////////////////////////////////////////////////////////////////
// timing Test behavior of blas gemv
//////////////////////////////////////////////////////////////////////////////
#include <nt2/toolbox/algebra/include/functions/gemv.hpp>
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

template<class T> struct dispatch_gemv_test
{
  typedef nt2::memory::
          container <nt2::tag::table_, nt2::id_<0>, T, nt2::settings(nt2::_2D)>  table_2D_type;
  typedef nt2::memory::
          container <nt2::tag::table_, nt2::id_<0>, T, nt2::settings(nt2::_1D)>  table_1D_type;

  dispatch_gemv_test( long int m, long int n, T alpha, T beta) 
    : a(of_size(m,n)), b(of_size(n)), r(of_size(n))
    , AL(alpha), BE(beta) 
  {
    dim1_a = size(a)(1);
    dim2_a = size(a)(2);
    
    for(std::size_t i = 1; i <= dim1_a; i++)
      for(std::size_t j = 1; j <= dim2_a; j++)
      {
        pos1 = boost::fusion::make_vector(i,j);
        a[pos1] =i*j;
      }

    for(std::size_t j = 1; j <= dim2_a; j++)
    {
      pos2 = boost::fusion::make_vector(j);
      b[pos2] = j;
    }

    for(std::size_t j = 1; j <= dim2_a; j++)
    {
      pos2 = boost::fusion::make_vector(j);
      r[pos2] = 0.0;
    }
  }
  
  void operator()()
  {
    nt2::gemv('N', a, b, r, AL, BE);
  }

  table_2D_type a;
  table_1D_type b, r;
  boost::fusion::vector<int,int> pos1;
  boost::fusion::vector<int> pos2;
  T AL, BE;
  std::size_t dim1_a;
  std::size_t dim2_a;
};

template<class T> struct system_gemv_test
{
  typedef nt2::memory::
          container <nt2::tag::table_, nt2::id_<0>, T, nt2::settings(nt2::_2D)>  table_2D_type;
  typedef nt2::memory::
          container <nt2::tag::table_, nt2::id_<0>, T, nt2::settings(nt2::_1D)>  table_1D_type;

  typedef typename table_2D_type::parent::lead_t lead_t;

  system_gemv_test(  const char c, long int m
                   , long int n, T alpha, T beta) 
    : a(of_size(m,n)), b(of_size(n)), r(of_size(n))
    , M(nt2::size(a)(c=='T'?2:1)), N(nt2::size(a)(c=='T'?1:2))
    , LDA(boost::simd::memory::align_on(size(a, 1), lead_t::value))
    , INCX(1)
    , INCY(1)
    , C(c), AL(alpha), BE(beta) 
  {
    dim1_a = size(a)(1);
    dim2_a = size(a)(2);
 
    for(std::size_t i = 1; i <= dim1_a; i++)
      for(std::size_t j = 1; j <= dim2_a; j++)
      {
        pos1 = boost::fusion::make_vector(i,j);
        a[pos1] =i*j;
      }

    for(std::size_t j = 1; j <= dim2_a; j++)
    {
      pos2 = boost::fusion::make_vector(j);
      b[pos2] = j;
    }

    for(std::size_t j = 1; j <= dim2_a; j++)
    {
      pos2 = boost::fusion::make_vector(j);
      r[pos2] = 0.0;
    }
  }
  
  void operator()()
  {
    nt2::ext::gemv(&C, &M, &N, &AL, a.begin(), &LDA, b.begin(), &INCX, &BE, r.begin(), &INCY);
  }

  table_2D_type a;
  table_1D_type b, r;
  long int M, N, K, LDA, INCX, INCY;
  boost::fusion::vector<int,int> pos1;
  boost::fusion::vector<int> pos2;
  T AL, BE;
  char C;
  std::size_t dim1_a;
  std::size_t dim2_a;
};

template<class T> struct normal_gemv_test
{
  typedef nt2::memory::
          container <nt2::tag::table_, nt2::id_<0>, T, nt2::settings(nt2::_2D)>  table_2D_type;
  typedef nt2::memory::
          container <nt2::tag::table_, nt2::id_<0>, T, nt2::settings(nt2::_1D)>  table_1D_type;
  typedef typename table_2D_type::parent::lead_t lead_t;

  normal_gemv_test(long int m, long int n) 
    : a(of_size(m,n)), b(of_size(n)), r(of_size(n))
  {
    dim1_a = size(a)(1);
    dim2_a = size(a)(2);
    
    for(std::size_t i = 1; i <= dim1_a; i++)
      for(std::size_t j = 1; j <= dim2_a; j++)
      {
        pos1 = boost::fusion::make_vector(i,j);
        a[pos1] =i*j;
      }

    for(std::size_t j = 1; j <= dim2_a; j++)
    {
      pos2 = boost::fusion::make_vector(j);
      b[pos2] = j;
    }

    for(std::size_t j = 1; j <= dim2_a; j++)
    {
      pos2 = boost::fusion::make_vector(j);
      r[pos2] = 0.0;
    }
  }
  
  void operator()()
  {
    for(std::size_t i = 1; i <= dim1_a; i++)
      for(std::size_t k = 1; k <= dim2_a; k++)
      {
        pos3 = boost::fusion::make_vector(i);
        pos2 = boost::fusion::make_vector(k);
        pos1 = boost::fusion::make_vector(i,k);
        r[pos3] += a[pos1]*b[pos2];
      } 
  }

  table_2D_type a; 
  table_1D_type b, r;
  boost::fusion::vector<int,int> pos1;
  boost::fusion::vector<int> pos2,pos3;
  std::size_t dim1_a;
  std::size_t dim2_a;
};

template<class T> void do_test()
{
  std::cout << "============================\n";
  std::cout << "Matrix product dummy\n";
  std::cout << "============================\n";
  normal_gemv_test<T> ngt(1024, 512);
  double dv1 = nt2::unit::perform_benchmark( ngt, 1.);
  std::cout << dv1/(1024*1024) << "\t\n";
 
  std::cout << "============================\n";
  std::cout << "Dispatch blas\n";
  std::cout << "============================\n";
  dispatch_gemv_test<T> dgt(1024, 512, 1.0, 0.0);
  double dv3 = nt2::unit::perform_benchmark( dgt, 1.);
  std::cout << dv3/(1024*1024) << "\t\n";

  std::cout << "============================\n";
  std::cout << "System blas\n";
  std::cout << "============================\n";
  system_gemv_test<T> sgt('N', 1024, 512, 1.0, 0.0);
  double dv2 = nt2::unit::perform_benchmark( sgt, 1.);
  std::cout << dv2/(1024*1024) << "\t\n";
}

NT2_TEST_CASE_TPL( gemv, (double)(float) )
{
  do_test<T>();
}
