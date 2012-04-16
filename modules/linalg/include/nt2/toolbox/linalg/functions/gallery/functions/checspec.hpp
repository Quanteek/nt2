/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_FUNCTIONS_CHEBSPEC_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_FUNCTIONS_CHEBSPEC_HPP_INCLUDED
#include <nt2/include/constants/pi.hpp>
#include <nt2/include/constants/zero.hpp>
#include <nt2/include/constants/one.hpp>
#include <nt2/include/constants/two.hpp>
#include <nt2/include/functions/ones.hpp>
#include <nt2/include/functions/zeros.hpp>
#include <nt2/include/functions/eye.hpp>
#include <nt2/include/functions/sqr.hpp>
#include <nt2/include/functions/oneminus.hpp>
#include <nt2/include/functions/colvect.hpp>
#include <nt2/include/functions/rowvect.hpp>
#include <nt2/include/functions/colon.hpp>
#include <nt2/include/functions/height.hpp>
#include <nt2/include/functions/width.hpp>
#include <nt2/table.hpp>


namespace nt2
{    //////////////////////////////////////////////////////////////////////////////
    //     chebspec  Chebyshev spectral differentiation matrix.
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( nt2::tag::magic_, tag::cpu_,
                                     (A0)(A1)(T), 
                                     (scalar_<integer_<A0> > )
                                     (scalar_<integer_<A1> > )
                                     ((target_<scalar_<floating_<T> > >)) 
                                     )
  {
    typedef T::type value_type; 
    typedef table<value_type> result_type; 
    inline result_type 
      operator()(size_t n, size_t k)
      {
        typedef typename table<value_type, C_index_> table_t;      
        if(k == 1) n++;
        n--;
        table_t c = zeros(n+1, n+1, T());
        table_t one = ones(n+1,1,T());
        table_t tone =  ones(1, n+1,T());
        value_type k1 = Pi<value_type>()/static_cast<value_type>(n); 
        table_t x = nt2::cos(nt2::colvect(nt2::colon(Zero<value_type>(),
                                              One<value_type>(),
                                              n))*k1);
        table_t d = ones(n+1,1, T());
        d(0) = Two<value_type>(); d(n) =  Two<value_type>();
        // eye(size(C)) on next line avoids div by zero.
        c = div( d*rowvect(rec(d)), (x*tone)-one*rowvect(x) + eye(size(c)));
        //  Now fix diagonal and signs.
        c(0) = (Two<value_type>()*sqr(n)+One<value_type>())/Six<value_type>();
        table<size_t, C_index_> rr, cc;
        rr =  ri(height(c), width(c))%2;
        cc =  ci(height(c), width(c))%2; 
        //      c = where(x_or(rr, cc), -c, c); 
        for (size_t i=1; i <n;  i++){
          c(i,i) = -x(i)/(Two<value_type>()*(oneminus(sqr(x(i)))));
        }
        c(n,n) = -c(0);
        if (k == 1){
          table_t c1 =  c(_(1, n), _(1, n));
          return c1; 
        }
        return c; 
      }
  };    

}


#endif

// /////////////////////////////////////////////////////////////////////////////
// End of chebspec.hpp
// /////////////////////////////////////////////////////////////////////////////
