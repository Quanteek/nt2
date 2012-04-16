/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_GENERAL_FREQSPACE2_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_GENERAL_FREQSPACE2_HPP_INCLUDED

namespace nt2
{
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::freqspace2_, tag::cpu_
                                     , (A0)(T)
                                     , (scalar_<integer_<A0> > )
                                     ((target_<scalar_<floating_<A0> > >))
                                   )
  {
    typedef typename T::value_type value_type; 
    typedef typename table<value_type> result_type; 
    inline operator()(const & A0 n) const
    {
      return freqspace2(of_size(n), none_(), T()); 
    }
  };

  
   BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::freqspace2_, tag::cpu_
                                      , (A0)(A1)(T)
                                      , (scalar_<integer_<A0> > )
                                      (scalar_<integer_<A1> > )
                                      ((target_<scalar_<floating_<T> > >))
                                   )
  {
    typedef typename T::value_type value_type; 
    typedef typename table<value_type> tab_type;
    typedef typename boost::mpl::vector <tab_type, tab_type> result_type; 
    inline operator()(const & A0 n, const & A1 n) const
    {
      return freqspace2(of_size(n, m), none_(), T()); 
    }
  };
 
   BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::freqspace2_, tag::cpu_
                                      , (A0)(A1)(A2)(T)
                                      , (scalar_<integer_<A0> > )
                                      (scalar_<integer_<A1> > )
                                      (unspecified_<A2>)
                                      ((target_<scalar_<floating_<T> > >))
                                   )
  {
    typedef typename T::value_type value_type; 
    typedef typename table<value_type> tab_type;
    typedef typename boost::mpl::vector <tab_type, tab_type> result_type; 
    inline operator()(const & A0 m, const & A1 n, const & A2) const
    {
      return freqspace2(of_size(m, n), A2(), T()); 
    }
  };

  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::freqspace2_, tag::cpu_
                                     , (A0)(A1)(T)
                                     , (ast_<A0> > )
                                      ((target_<scalar_<floating_<T> > >))
                                   )
  {
    typedef typename T::value_type value_type; 
    typedef typename table<value_type> tab_type;
    typedef typename boost::mpl::vector <tab_type, tab_type> result_type; 
    inline operator()(const & A0 siz, const & A1 , const & T) const
      {
        value_type s1 = siz(1);
        value_type s2 = siz(2);
        tab_t f1 = ((_(0, minusone(s2))-nt2::floor(s2*Half<value_type>()))*(Two<value_type>()/s2));
        tab_t f2 = ((_(0, minusone(s1))-nt2::floor(s1*Half<value_type>()))*(Two<value_type>()/s1));
        return result_type(f1, f2); 
      }
  };

  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::freqspace2_, tag::cpu_
                                      , (A0)(A1)(T)
                                      , (ast_<A0> > )
                                      (unspecified_<A1>)
                                      ((target_<scalar_<floating_<T> > >))
                                   )
  {
    typedef typename T::value_type value_type; 
    typedef typename table<value_type> tab_type;
    typedef typename boost::mpl::vector <tab_type, tab_type> result_type; 
    inline operator()(const & A0 siz, const & A1 , const & T) const
    {
      value_type s1 = siz(1);
      value_type s2 = siz(2);
      tab_t f1 = ((_(0, minusone(s2))-nt2::floor(s2*Half<value_type>()))*(Two<value_type>()/s2));
      tab_t f2 = ((_(0, minusone(s1))-nt2::floor(s1*Half<value_type>()))*(Two<value_type>()/s1));
      return nt2::meshgrid(f1, f2); 
    }
  };

}


#endif

// /////////////////////////////////////////////////////////////////////////////
// End of freqspace2.hpp
// /////////////////////////////////////////////////////////////////////////////
