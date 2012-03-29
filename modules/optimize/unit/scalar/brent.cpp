//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II       
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI       
//                                                                              
//          Distributed under the Boost Software License, Version 1.0.          
//                 See accompanying file LICENSE.txt or copy at                 
//                     http://www.boost.org/LICENSE_1_0.txt                     
//==============================================================================
#define NT2_UNIT_MODULE "nt2 optimize toolbox - brent"

//////////////////////////////////////////////////////////////////////////////
// unit test behavior of optimize components
//////////////////////////////////////////////////////////////////////////////
#include <nt2/include/functions/brent.hpp>
#include <nt2/sdk/unit/tests.hpp>
#include <nt2/sdk/unit/module.hpp>
#include <iostream>

  template < class T > 
  T f0(const T& x, T a,  T b) {
    return x*x*x-3*x+4;
  }
  
  template < class T > 
  void f00(T & y, const T& x, T a,  T b)      {
    y = x*x*x-a*x+b; 
  }
  
  template < class T > 
  class fp : public nt2::optimizableProc < T > 
  {
  public:  
    fp(T a_, T b_): a(a_), b(b_){}; 
    virtual void operator()(T& y, const T& x ) const
    {
      return f00(y, x, a, b); 
    }
  private:
    T a;
    T b; 
  };
  
  template < class T > 
  class f : public nt2::optimizableFunc < T >
  {
  public:  
    f(T a_, T b_): a(a_), b(b_){};  
    virtual T operator()( const T& x ) const
    {
      return f0(x, a, b); 
    }
  private:
    T a;
    T b; 
  };
  
  double f1d(const double& x )     {
    return x*x*x-3*x+4;
  }
  
  void f2d(double & y, const double& x )   { 
    y = x*x*x-3*x+4;
  }

  float f1f(const float& x )     {
    return x*x*x-3*x+4;
  }
  
  void f2f(float & y, const float& x )   { 
    y = x*x*x-3*x+4;
  }

NT2_TEST_CASE_TPL ( brentd, (double)) 
{
  using nt2::brent;
  
  brent < T > br;
  fp<T> proc2minimize(3, 4); 
  br.optimize(proc2minimize, T(0), T(0.5), T(2));
  
  std::cout << "minimum " << br.getMinimum() <<  " au point " <<  br.getMinimumPosition()
            << " en " << br.getNbIteration() <<  " iterations" << std::endl;
  NT2_TEST(nt2::abs(br.getMinimumPosition()- T(1)) <= nt2::Sqrteps<T>());
  
  br.optimize(&f2d, T(0), T(0.5), T(2)); 
  std::cout << "minimum " << br.getMinimum() <<  " au point " <<  br.getMinimumPosition()
            << " en " << br.getNbIteration() <<  " iterations" << std::endl; 
  NT2_TEST(nt2::abs(br.getMinimumPosition()- T(1)) <= nt2::Sqrteps<T>());
  
  f<T> func2minimize(3, 4); 
  br.optimize(func2minimize, T(0), T(0.5), T(2)); 
  std::cout << "minimum " << br.getMinimum() <<  " au point " <<  br.getMinimumPosition()
            << " en " << br.getNbIteration() <<  " iterations" << std::endl; 
  NT2_TEST(nt2::abs(br.getMinimumPosition()- T(1)) <= nt2::Sqrteps<T>());
  
  br.optimize(&f1d, T(0), T(0.5), T(2)); 
  std::cout << "minimum " << br.getMinimum() <<  " au point " <<  br.getMinimumPosition()
            << " en " << br.getNbIteration() <<  " iterations" << std::endl; 
  NT2_TEST(nt2::abs(br.getMinimumPosition()- T(1)) <= nt2::Sqrteps<T>());
}
NT2_TEST_CASE_TPL ( brentf, (float)) 
{
  using nt2::brent;
  
  brent < T > br;
  fp<T> proc2minimize(3, 4); 
  br.optimize(proc2minimize, T(0), T(0.5), T(2)); 
  std::cout << "minimum " << br.getMinimum() <<  " au point " <<  br.getMinimumPosition()
            << " en " << br.getNbIteration() <<  " iterations" << std::endl; 
  NT2_TEST(nt2::abs(br.getMinimumPosition()- T(1)) <= nt2::Sqrteps<T>());
  
  br.optimize(&f2f, T(0), T(0.5), T(2)); 
  std::cout << "minimum " << br.getMinimum() <<  " au point " <<  br.getMinimumPosition()
            << " en " << br.getNbIteration() <<  " iterations" << std::endl; 
  NT2_TEST(nt2::abs(br.getMinimumPosition()- T(1)) <= nt2::Sqrteps<T>());
  
  f<T> func2minimize(3, 4); 
  br.optimize(func2minimize, T(0), T(0.5), T(2)); 
  std::cout << "minimum " << br.getMinimum() <<  " au point " <<  br.getMinimumPosition()
            << " en " << br.getNbIteration() <<  " iterations" << std::endl; 
  NT2_TEST(nt2::abs(br.getMinimumPosition()- T(1)) <= nt2::Sqrteps<T>());
  
  br.optimize(&f1f, T(0), T(0.5), T(2)); 
  std::cout << "minimum " << br.getMinimum() <<  " au point " <<  br.getMinimumPosition()
            << " en " << br.getNbIteration() <<  " iterations" << std::endl; 
  NT2_TEST(nt2::abs(br.getMinimumPosition()- T(1)) <= nt2::Sqrteps<T>());
}
