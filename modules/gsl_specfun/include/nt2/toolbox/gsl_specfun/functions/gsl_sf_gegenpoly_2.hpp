//////////////////////////////////////////////////////////////////////////////
///   Copyright 2003 and onward LASMEA UMR 6602 CNRS/U.B.P Clermont-Ferrand
///   Copyright 2009 and onward LRI    UMR 8623 CNRS/Univ Paris Sud XI
///
///          Distributed under the Boost Software License, Version 1.0
///                 See accompanying file LICENSE.txt or copy at
///                     http://www.boost.org/LICENSE_1_0.txt
//////////////////////////////////////////////////////////////////////////////
/*!
 * \file
**/
#ifndef NT2_TOOLBOX_GSL_SPECFUN_FUNCTIONS_GSL_SF_GEGENPOLY_2_HPP_INCLUDED
#define NT2_TOOLBOX_GSL_SPECFUN_FUNCTIONS_GSL_SF_GEGENPOLY_2_HPP_INCLUDED
#include <nt2/include/simd.hpp>
#include <nt2/include/functor.hpp>

/*!
 * \ingroup gsl_specfun
 * \defgroup gsl_specfun_gsl_sf_gegenpoly_2 gsl_sf_gegenpoly_2
 *
 * \par Description
 * TODO Put description here
 *
 * \par Header file
 * 
 * \code
 * #include <nt2/toolbox/gsl_specfun/include/functions/gsl_sf_gegenpoly_2.hpp>
 * \endcode
 * 
 * 
 * \synopsis
 *
 * \code
 * namespace nt2
 * {
 *   namespace gsl_specfun
 *   {
 *     template <class A0>
 *       meta::call<tag::gsl_sf_gegenpoly_2_(A0,A0)>::type
 *       gsl_sf_gegenpoly_2(const A0 & a0,const A0 & a1);
 *   }
 * }
 * \endcode
 *
 * \param a0 the first parameter of gsl_sf_gegenpoly_2
 * \param a1 the second parameter of gsl_sf_gegenpoly_2
 * 
 * \return a value of the common type of the parameters
 *  
 * \par Notes
 * In SIMD mode, this function acts elementwise on the inputs vectors elements
 * \par
 * \par
 * gsl_specfun library defines functions for float and double entries.
 * \par
 * As they are written in C the original name of the float version is
 * generally terminated by and extra 'f',
 * this is not the case for the nt2 version which dispatch to
 * the correct function according to the inputs types.
 * \par
 * gsl_specfun library defines functions for double entries only.
 * Nevertheless, they can be called with float entries under nt2 calls
 * to return float outputs.
 *  
**/

namespace nt2 { namespace gsl_specfun { namespace tag
  {         
    /*!
     * \brief Define the tag gsl_sf_gegenpoly_2_ of functor gsl_sf_gegenpoly_2 
     *        in namespace nt2::gsl_specfun::tag for toolbox gsl_specfun
    **/
    struct gsl_sf_gegenpoly_2_ : ext::elementwise_<gsl_sf_gegenpoly_2_> { typedef ext::elementwise_<gsl_sf_gegenpoly_2_> parent; };
  }
  NT2_FUNCTION_IMPLEMENTATION(gsl_specfun::tag::gsl_sf_gegenpoly_2_, gsl_sf_gegenpoly_2, 2)
  } }

#include <nt2/toolbox/gsl_specfun/functions/scalar/gsl_sf_gegenpoly_2.hpp>
// #include <nt2/toolbox/gsl_specfun/functions/simd/all/gsl_sf_gegenpoly_2.hpp> 

#endif

// modified by jt the 29/12/2010
