/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_GENERAL_FREQSPACE1_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_GENERAL_FREQSPACE1_HPP_INCLUDED

namespace nt2
{
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::freqspace1_, tag::cpu_
                                     , (A0)(T)
                                     , (scalar_<integer_<A0> > )
                                     ((target_<scalar_<floating_<A0> > >))
                                   )
  {
    typedef typename T::value_type value_type; 
    typedef typename table<value_type> result_type; 
    inline operator()(const & A0 n) const
    {
      return  n == 0
        ? zeros(1, 0, T())
        :_(Zero<value_type>(), Two<value_type>/n:One<value_type>(), One<value_type>());
    }
  };

  
   BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::freqspace1_, tag::cpu_
                                      , (A0)(A1)(T)
                                      , (scalar_<integer_<A0> > )
                                      (unspecified_<A1>)
                                      ((target_<scalar_<floating_<T> > >))
                                   )
  {
    typedef typename T::value_type value_type; 
    typedef typename table<value_type> result_type; 
    inline operator()(const & A0 m, const & A1) const
    {
      value_type tmp = Two<value_type>/n; 
      return n == 0
        ? zeros(1, 0, T())
        :_(Zero<value_type>(), tmp, tmp*(n-1));
    }
  };


    

}


#endif
