/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_DETAILS_UTILITY_LOWER_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_DETAILS_UTILITY_LOWER_HPP_INCLUDED



namespace nt2
{
  namespace details
  {
    
    inline char lower(const char c){
      return (c >= 'A' && c <= 'Z') ? c +('a'-'A') : c; 
    }
    inline char upper(const char c){
      return (c >= 'a' && c <= 'z') ? c -('a'-'A') : c; 
    }
  }
}


#endif

// /////////////////////////////////////////////////////////////////////////////
// End of lower.hpp
// /////////////////////////////////////////////////////////////////////////////
