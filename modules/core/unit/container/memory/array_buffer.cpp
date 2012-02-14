//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#define NT2_UNIT_MODULE "nt2::memory boost::array adaptation as a Buffer"

#include <nt2/sdk/memory/array_buffer.hpp>

#include <algorithm>
#include <boost/array.hpp>
#include <boost/fusion/adapted/array.hpp>
#include <boost/fusion/include/vector_tie.hpp>
#include <boost/fusion/include/make_vector.hpp>

#include <nt2/sdk/unit/module.hpp>
#include <nt2/sdk/unit/tests/basic.hpp>
#include <nt2/sdk/unit/tests/relation.hpp>
#include <nt2/sdk/unit/tests/type_expr.hpp>
#include <nt2/sdk/unit/tests/exceptions.hpp>

//==============================================================================
// array_buffer value_of traits
//==============================================================================
NT2_TEST_CASE_TPL( array_buffer_value_of, NT2_TYPES )
{
  using boost::dispatch::meta::value_of;
  using nt2::memory::array_buffer;
  using boost::mpl::_;

  typedef array_buffer<T,2,0> buffer_type;

  NT2_TEST_EXPR_TYPE(buffer_type(), value_of<_>, T );
}

//==============================================================================
// array_buffer model_of traits
//==============================================================================
template<class T, class U>
struct  apply_model
      : boost::mpl::apply < typename boost::dispatch::meta::model_of<T>::type
                          , U
                          >
{};

NT2_TEST_CASE_TPL( array_buffer_model_of, NT2_TYPES )
{
  using boost::dispatch::meta::model_of;
  using nt2::memory::array_buffer;
  using boost::mpl::_;

  typedef array_buffer<char,2,0> buffer_type;

  NT2_TEST_EXPR_TYPE((buffer_type()),(apply_model<_,T>),(array_buffer<T,2,0>));
}

//==============================================================================
// Test for dynamic default array_buffer ctor
//==============================================================================
NT2_TEST_CASE_TPL( array_buffer_default_ctor, NT2_TYPES)
{
  using nt2::memory::array_buffer;

  array_buffer<T,2,1> b;

  NT2_TEST_EQUAL(b.size()       , 2 );
  NT2_TEST_EQUAL(b.inner_size() , 2 );
  NT2_TEST_EQUAL(b.outer_size() , 1 );
  NT2_TEST_EQUAL(b.lower()      , 1 );
  NT2_TEST_EQUAL(b.inner_lower(), 1 );
  NT2_TEST_EQUAL(b.outer_lower(), 1 );
  NT2_TEST_EQUAL(b.upper()      , 2 );
  NT2_TEST_EQUAL(b.inner_upper(), 2 );
  NT2_TEST_EQUAL(b.outer_upper(), 1 );
}

//==============================================================================
// Test for dynamic array_buffer copy ctor
//==============================================================================
NT2_TEST_CASE_TPL( array_buffer_data_ctor, NT2_TYPES)
{
  using nt2::memory::array_buffer;

  typedef array_buffer<T,5,-2> buffer_type ;
  buffer_type b;

  for ( typename buffer_type::difference_type i = b.lower(); i <= b.upper(); ++i )
    b(i) = T(3+i);

  buffer_type x(b);

  NT2_TEST_EQUAL(x.size()       ,  5 );
  NT2_TEST_EQUAL(x.inner_size() ,  5 );
  NT2_TEST_EQUAL(x.outer_size() ,  1 );
  NT2_TEST_EQUAL(x.lower()      , -2 );
  NT2_TEST_EQUAL(x.inner_lower(), -2 );
  NT2_TEST_EQUAL(x.outer_lower(),  1 );
  NT2_TEST_EQUAL(x.upper()      ,  2 );
  NT2_TEST_EQUAL(x.inner_upper(),  2 );
  NT2_TEST_EQUAL(x.outer_upper(),  1 );

  for ( typename buffer_type::difference_type i = x.lower(); i <= x.upper(); ++i )
    NT2_TEST_EQUAL( x(i), T(3+i) );
}

//==============================================================================
// Test for array_buffer asserting constructor
//==============================================================================
NT2_TEST_CASE_TPL( array_buffer_wrong_data_ctor, NT2_TYPES)
{
  using nt2::memory::array_buffer;
  typedef array_buffer<T,5,-2> buffer_type ;

  NT2_TEST_ASSERT( buffer_type too_much(7) );
  NT2_TEST_ASSERT( buffer_type too_few(2)  );
}

//==============================================================================
// Test for array_buffer assignment
//==============================================================================
NT2_TEST_CASE_TPL(array_buffer_assignment, NT2_TYPES )
{
  using nt2::memory::array_buffer;

  typedef array_buffer<T,5,-2> buffer_type ;
  buffer_type x, b;

  for ( typename buffer_type::difference_type i = b.lower(); i <= b.upper(); ++i )
    b(i) = T(3+i);

  x = b;

  NT2_TEST_EQUAL(x.size()       ,  5 );
  NT2_TEST_EQUAL(x.inner_size() ,  5 );
  NT2_TEST_EQUAL(x.outer_size() ,  1 );
  NT2_TEST_EQUAL(x.lower()      , -2 );
  NT2_TEST_EQUAL(x.inner_lower(), -2 );
  NT2_TEST_EQUAL(x.outer_lower(),  1 );
  NT2_TEST_EQUAL(x.upper()      ,  2 );
  NT2_TEST_EQUAL(x.inner_upper(),  2 );
  NT2_TEST_EQUAL(x.outer_upper(),  1 );

  for ( typename buffer_type::difference_type i = x.lower(); i <= x.upper(); ++i )
    NT2_TEST_EQUAL( x(i), T(3+i) );
}

//==============================================================================
// Test for array_buffer swap
//==============================================================================
NT2_TEST_CASE_TPL(array_buffer_swap, NT2_TYPES )
{
  using nt2::memory::array_buffer;

  typedef array_buffer<T,5,1> buffer_type ;

  buffer_type b;
  buffer_type x;

  for ( typename buffer_type::difference_type i = b.lower(); i <= b.upper(); ++i )
    b(i) = T(3+i);

  for ( typename buffer_type::difference_type i = x.lower(); i <= x.upper(); ++i )
    x(i) = T(10*i);

  swap(b,x);

  NT2_TEST_EQUAL(x.size()       ,  5 );
  NT2_TEST_EQUAL(x.inner_size() ,  5 );
  NT2_TEST_EQUAL(x.outer_size() ,  1 );
  NT2_TEST_EQUAL(x.lower()      ,  1 );
  NT2_TEST_EQUAL(x.inner_lower(),  1 );
  NT2_TEST_EQUAL(x.outer_lower(),  1 );
  NT2_TEST_EQUAL(x.upper()      ,  5 );
  NT2_TEST_EQUAL(x.inner_upper(),  5 );
  NT2_TEST_EQUAL(x.outer_upper(),  1 );

  NT2_TEST_EQUAL(b.size()       ,  5 );
  NT2_TEST_EQUAL(b.inner_size() ,  5 );
  NT2_TEST_EQUAL(b.outer_size() ,  1 );
  NT2_TEST_EQUAL(b.lower()      ,  1 );
  NT2_TEST_EQUAL(b.inner_lower(),  1 );
  NT2_TEST_EQUAL(b.outer_lower(),  1 );
  NT2_TEST_EQUAL(b.upper()      ,  5 );
  NT2_TEST_EQUAL(b.inner_upper(),  5 );
  NT2_TEST_EQUAL(b.outer_upper(),  1 );

  for ( typename buffer_type::difference_type i = x.lower(); i <= x.upper(); ++i )
    NT2_TEST_EQUAL( x(i), T(3+i) );

  for ( typename buffer_type::difference_type i = b.lower(); i <= b.upper(); ++i )
    NT2_TEST_EQUAL( b(i), T(10*i) );
}

//==============================================================================
// array_buffer Range interface
//==============================================================================
struct f_
{
  template<class T> T operator()(T const& e) const { return T(10*e); }
};

NT2_TEST_CASE_TPL(array_buffer_iterator, NT2_TYPES )
{
  using nt2::memory::array_buffer;

  typedef array_buffer<T,5,-2> buffer_type ;

  buffer_type x;
  for ( typename buffer_type::difference_type i = x.lower(); i <= x.upper(); ++i )
  x(i) = T(3+i);

  f_ f;

  typename buffer_type::iterator b = x.begin();
  typename buffer_type::iterator e = x.end();

  std::transform(b,e,b,f);

  for ( typename buffer_type::difference_type i = x.lower(); i <= x.upper(); ++i )
    NT2_TEST_EQUAL( x(i), T(f(3+i)) );
}

//==============================================================================
// Test for array_buffer asserting resize
//==============================================================================
NT2_TEST_CASE_TPL( array_buffer_wrong_resize, NT2_TYPES)
{
  using nt2::memory::array_buffer;
  typedef array_buffer<T,5,-2> buffer_type ;
  buffer_type b;

  NT2_TEST_ASSERT( b.resize(7) );
  NT2_TEST_ASSERT( b.resize(1) );
}
