/*******************************************************************************
 *         Copyright 2003 & onward LASMEA UMR 6602 CNRS/Univ. Clermont II
 *         Copyright 2009 & onward LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#define NT2_UNIT_MODULE "nt2 container container"

#include <nt2/sdk/memory/container.hpp>
#include <nt2/core/functions/of_size.hpp>
#include <nt2/core/container/table/normalize_settings.hpp>
#include <boost/fusion/include/vector_tie.hpp>

#include <iostream>
#include <nt2/sdk/timing/now.hpp>
#include <nt2/sdk/unit/details/helpers.hpp>
#include <nt2/sdk/unit/perform_benchmark.hpp>
#include <nt2/sdk/unit/module.hpp>

template<class T> struct container_1D_dynamic_test
{
  typedef nt2::memory::
          container <nt2::tag::table_, nt2::id_<0>, T, nt2::settings(nt2::_1D)>  buffer_t;

  container_1D_dynamic_test ( std::ptrdiff_t s0 )
                            : data(boost::fusion::vector_tie(s0))
                            , data2(boost::fusion::vector_tie(s0))
                            , s0_(s0)
  {}

  void operator()()
  {
    for(std::ptrdiff_t i = 1; i <= s0_; ++i)
      data[boost::fusion::vector_tie(i)] = data2[boost::fusion::vector_tie(i)];
  }

  buffer_t        data,data2;
  std::ptrdiff_t  s0_;
};

template<class T> struct container_1D_static_test
{
  typedef nt2::memory::
          container < nt2::tag::table_, nt2::id_<0>
                    , T
                    , nt2::settings(nt2::of_size_<256*256>, nt2::automatic_)
                    >  buffer_t;


  container_1D_static_test( std::ptrdiff_t s0 )
                          : data(boost::fusion::vector_tie(s0))
                          , data2(boost::fusion::vector_tie(s0))
                          , s0_(s0)
  {}

  void operator()()
  {
    for(std::ptrdiff_t i = 1; i <= s0_; ++i)
      data[boost::fusion::vector_tie(i)] = data2[boost::fusion::vector_tie(i)];
  }

  buffer_t        data,data2;
  std::ptrdiff_t  s0_;
};

template<class T> struct container_2D_static_test
{
 typedef nt2::memory::
          container < nt2::tag::table_, nt2::id_<0>, T
                    , nt2::settings(nt2::of_size_<256,256>, nt2::automatic_)
                    >  buffer_t;

  container_2D_static_test ( std::size_t s0, std::size_t s1 )
                              : data( boost::fusion::vector_tie(s0,s1) )
                              , data2( boost::fusion::vector_tie(s0,s1) )
                              , s0_(s0), s1_(s1)
  {}

  void operator()()
  {
    for(std::ptrdiff_t j = 1; j <= s1_; ++j)
      for(std::ptrdiff_t i = 1; i <= s0_; ++i)
        data[boost::fusion::vector_tie(i,j)] = data2[boost::fusion::vector_tie(i,j)];
  }

  buffer_t        data,data2;
  std::ptrdiff_t  s0_,s1_;
};


template<class T> struct container_2D_dynamic_test
{
  typedef nt2::memory::
          container <nt2::tag::table_, nt2::id_<0>, T, nt2::settings()>  buffer_t;

  container_2D_dynamic_test ( std::size_t s0, std::size_t s1 )
                              : data( boost::fusion::vector_tie(s0,s1) )
                              , data2( boost::fusion::vector_tie(s0,s1) )
                              , s0_(s0), s1_(s1)
  {}

  void operator()()
  {
    for(std::ptrdiff_t j = 1; j <= s1_; ++j)
      for(std::ptrdiff_t i = 1; i <= s0_; ++i)
        data[boost::fusion::vector_tie(i,j)] = data2[boost::fusion::vector_tie(i,j)];
  }

  buffer_t        data,data2;
  std::ptrdiff_t  s0_,s1_;
};

template<class T> struct std_1D_test
{
  std_1D_test(std::ptrdiff_t h, std::ptrdiff_t w) : data(h*w),data2(h*w)
  {}

  void operator()()
  {
    for(std::size_t i = 0; i < data.size(); ++i)
      data[i] = data2[i];
  }

  std::vector<T,boost::simd::memory::allocator<T> > data,data2;
};

template<class T> struct std_2D_test
{
  std_2D_test(std::ptrdiff_t h, std::ptrdiff_t w) : data(h*w),data2(h*w)
  {
    s = h;
    d = w;
  }

  void operator()()
  {
    for(std::ptrdiff_t j = 0; j < d; ++j)
      for(std::ptrdiff_t i = 0; i < s; ++i)
        data[i+s*j] = data2[i+s*j];
  }

  std::ptrdiff_t d,s;
  std::vector<T,boost::simd::memory::allocator<T> > data,data2;
};

template<class T> void do_large(int H, int W)
{
  container_1D_dynamic_test<T> b(H*W);
  double d = nt2::unit::perform_benchmark( b, 1.)/2.;

  container_2D_dynamic_test<T> c(H,W);
  double e = nt2::unit::perform_benchmark( c, 1.)/2.;

  std_1D_test<T> z(H,W);
  double w = nt2::unit::perform_benchmark( z, 1.)/2.;

  std_2D_test<T> y(H,W);
  double v = nt2::unit::perform_benchmark( y, 1.)/2.;

  printf( "%d x %d : 1D %3.3f %3.3f (%3.3f%%) | 2D: %3.3f %3.3f (%3.3f%%)\n"
        , H, W
        , d/(H*W), w/(H*W), ((d-w)/w)*100
        , e/(H*W), v/(H*W), ((e-v)/v)*100
        );
}

NT2_TEST_CASE_TPL( container_large, (double)(float)(short)(char) )
{
  do_large<T>(320 , 240);
  do_large<T>(640 , 480);
  do_large<T>(1024, 768);
  do_large<T>(2048, 1556);
  do_large<T>(4096, 4096);
}

template<class T> void do_small(int H, int W)
{
  container_1D_dynamic_test<T> b(H*W);
  double d = nt2::unit::perform_benchmark( b, 1.)/2.;

  std_1D_test<T> z(H,W);
  double w = nt2::unit::perform_benchmark( z, 1.)/2.;

  container_2D_dynamic_test<T> c(H,W);
  double e = nt2::unit::perform_benchmark( c, 1.)/2.;

  container_1D_static_test<T> bs(H*W);
  double ds = nt2::unit::perform_benchmark( bs, 1.)/2.;

  container_2D_static_test<T> cs(H,W);
  double es = nt2::unit::perform_benchmark( cs, 1.)/2.;

  std_2D_test<T> y(H,W);
  double v = nt2::unit::perform_benchmark( y, 1.)/2.;

  std::cout << H << "x" << W << " : "
            << "1D "  << d/(H*W) << " (" <<  ((d-w)/w)*100    << ") "
                      << ds/(H*W) << " [" <<  ((ds-w)/w)*100  << "] "
                      << w/(H*W)
            << " | "
            << "2D "  << e/(H*W) << " (" <<  ((e-v)/v)*100    << ") "
                      << es/(H*W) << " [" <<  ((es-v)/v)*100  << "] "
                      << v/(H*W)
            << "\n";
}

NT2_TEST_CASE_TPL( container_small, (double)(float)(short)(char) )
{
  for(int N=1;N<=256;N*=2)
    do_small<T>(N,N);
}
