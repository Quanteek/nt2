/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_ALGEBRA_FUNCTIONS_PROPER_TRANS_HPP_INCLUDED
#define NT2_TOOLBOX_ALGEBRA_FUNCTIONS_PROPER_TRANS_HPP_INCLUDED

namespace nt2
{
  namespace ext
  {
    
    NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::trans_, tag::cpu_, 
                                (A)(M)(N)(SHA)(STA), 
                                ((table_< floating_<A>, nt2::settings(nt2::of_size_<N::value, M::value>,SHA,STA) > ))
                                )
    {
      typedef typename A::value_type value_type; 
      typedef nt2::table < value_type, nt2::settings(nt2:of_size_<N::value, M::value>, SHA, STA)> result_type;
      
      BOOST_FORCEINLINE result_type operator()(A const& a)
      {
        table < value_type > ta(N::value, M::value);
         
         for(std::size_t i=1; i <= size(a, 1); i++)
           {
             for(std::size_t j=1; j <= size(a, 2); j++)
               {
                 ta(i, j) = a(j, i); 
               }
           }
      }
    };
  }
}


#endif

// /////////////////////////////////////////////////////////////////////////////
// End of trans.hpp<2>
// /////////////////////////////////////////////////////////////////////////////
