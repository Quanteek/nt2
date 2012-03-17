/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_PASCAL_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_GALLERY_PASCAL_HPP_INCLUDED

namespace nt2
{
  template < class T > struct elmat
    {
      typedef table<T>         tab_t;

      static M_t pascal(const size_t &n,  const size_t & k = 0)
      {
        tab_t p(diag(2*(colon(1, n)%2)-1));
        p(_, 0) = 1;
        // Generate pascal(n,1): the Pascal Cholesky factor (up to signs).
        for(size_t j=1; j < n-1; j++)
          {

            for(size_t i=j+1; i < n; i++)
              {
                //      Range r(j+1, n-1);
                //      Range r1(j, n-2);
                p(i,j) = p(i-1,j) - p(i-1,j-1);
              }
          }
        if (k == 0 || k > 2)
          return prodMtM(p, p);
        else if (k == 2)
          return rot90(p, 3);
        else
          return p;
      }
    };

  template < class T >  inline
  typename pascal<T>::tab_t
  pascal(const size_t &n,  const size_t & k = 0)
  {
    return elmat<T>::pascal(n, k); 
  }
}


#endif

// /////////////////////////////////////////////////////////////////////////////
// End of pascal.hpp
// /////////////////////////////////////////////////////////////////////////////
