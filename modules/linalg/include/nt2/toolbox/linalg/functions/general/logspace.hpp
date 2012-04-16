/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_GENERAL_LOGSPACE_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_GENERAL_LOGSPACE_HPP_INCLUDED
// this matlab function is really bad design!
// the special pi case is quite absurd
namespace nt2
{
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::logspace_, tag::cpu_
                                     , (A0)(A1)
                                     , (scalar_<floating_<A0> > )
                                       (scalar_<floating_<A0> > )
                                       (scalar_<integer_<A1> > )
                                   )
  {
    typedef typename A0 value_type; 
    typedef typename table<value_type> result_type; 
    inline operator()(const & A0 d1, const & A0 dd2, const & A1 n) const
    {
      value_type d2 =  (dd2 == Pi<value_type>())
        ? dd2 = nt2::log10(Pi<value_type>());  //TODO define a constant
        : dd2;
      value_type fn =  value_type(n); 
      return  n
        ? exp10(cath(d1+_(Zero<value_type>(), fn-Two<value_type>())*(d2-d1)/(minusone(fn)), d2))
        : zeros(0, 1, as_<value_type>());
    }
  };

  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::logspace_, tag::cpu_
                                     , (A0)(A1)
                                     , (scalar_<floating_<A0> > )
                                       (scalar_<floating_<A0> > )
                                   )
  {
    typedef typename A0 value_type; 
    typedef typename table<value_type> result_type; 
    inline operator()(const & A0 d1, const & A0 dd2) const
    {
      return logspace(d1, d2, 50); 
    }
  };
 
end

    

}


#endif
