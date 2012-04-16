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
    
//   function [f1,f2] = freqspace2(n,flag)
// %FREQSPACE2 Frequency spacing for frequency response.
// %   FREQSPACE2 returns the implied frequency range for equally spaced
// %   frequency responses.  FREQSPACE2 is useful when creating desired
// %   frequency responses for FSAMP2, FWIND1, and FWIND2 as well as
// %   for various 1-D applications.
// %   
// %   [F1,F2] = FREQSPACE2(N) returns the 2-D frequency range vectors
// %   F1 and F2 for an N-by-N matrix.
// %   [F1,F2] = FREQSPACE2([M N]) returns the 2-D frequency range 
// %   vectors for an M-by-N matrix.
// %
// %   For 2-D vectors and n odd,  F = (-1+1/n:2/n:1-1/n).
// %   For 2-D vectors and n even, F = (-1    :2/n:1-2/n).
// %
// %   [F1,F2] = FREQSPACE2(...,'meshgrid') is equivalent to
// %       [F1,F2] = freqspace2(...); [F1,F2] = meshgrid(F1,F2);
// %
// %   F = FREQSPACE2(N) returns the 1-D frequency vector F assuming N
// %   equally spaced points around the unit circle.  For 1-D vectors, 
// %   F = (0:2/N:1).  F = FREQSPACE2(N,'whole') returns all N equally
// %   spaced points. In this case, F = (0:2/N:2*(N-1)/N).
// %
// %   Class support for inputs M,N:
// %      float: double, single
// %
// %   See also FSAMP2, FWIND1, FWIND2.
// %
// %   Note: FSAMP2, FWIND1 and FWIND2 are in the Image Processing Toolbox.

// %   Copyright 1984-2010 The MathWorks, Inc. 
// %   $Revision: 1.23.4.4 $  $Date: 2010/11/22 02:45:59 $

// if length(n)==1 && nargout>1, n = [n n]; end
// if nargin>1,
//   if ~ischar(flag),
//     error(message('MATLAB:freqspace2:Arg2NotStr'));
//   end
// end

// if nargout>1,
//   f1 = ((0:n(2)-1)-floor(n(2)/2))*(2/(n(2)));
//   f2 = ((0:n(1)-1)-floor(n(1)/2))*(2/(n(1)));
//   if nargin>1,
//     [f1,f2] = meshgrid(f1,f2);
//   end
// else
//   if nargin>1,
//     f1 = (0:2/n:2*(n-1)/n);
//   else
//     if length(n)==1 && n==0,
//       f1 = zeros(1,0,class(n));
//     else
//       f1 = (0:2/n:1);
//     end
//   end
// end
}


#endif

// /////////////////////////////////////////////////////////////////////////////
// End of freqspace2.hpp
// /////////////////////////////////////////////////////////////////////////////
