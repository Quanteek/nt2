//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II       
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI       
//                                                                              
//          Distributed under the Boost Software License, Version 1.0.          
//                 See accompanying file LICENSE.txt or copy at                 
//                     http://www.boost.org/LICENSE_1_0.txt                     
//==============================================================================
#ifndef NT2_TOOLBOX_LINALG_DETAILS_UTILITY_F77_WRAPPER_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_DETAILS_UTILITY_F77_WRAPPER_HPP_INCLUDED
////////////////////////////////////////////////////////////////////////////////
// This macro can be ifdefed to produce the correct symbol output for the
// FORTRAN library in use, adding or not the underscore to FORTRAN routine names
////////////////////////////////////////////////////////////////////////////////

#define NT2_F77NAME(x) BOOST_PP_CAT(x,_)
typedef long int la_int;
typedef void la_complex; 

#endif
