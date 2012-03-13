/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_DETAILS_LAPACK_UTILITY_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_DETAILS_LAPACK_UTILITY_HPP_INCLUDED
#include <string>
#include <nt2/toolbox/linalg/details/utility/f77_wrapper.hpp>
namespace nt2
{
  namespace details
  {
    extern "C"
    {
      int ilaenv_(const int *i, const char *n, const char *opts, const int *n1,
                  const int *n2, const int *n3, const int *n4,
                  long int n_len, long int opts_len);
    }
    ////////////////////////////////////////////////////////////////////////////
    // utility
    ////////////////////////////////////////////////////////////////////////////
    inline int EnvBlockSize(int ispec,
                            const std::string & fname,
                            const std::string & opts,
                            int n1 = -1, 
                            int n2 = -1,
                            int n3 = -1,
                            int n4 = -1)
    {
      int i = ispec;
      int N1 = n1;
      int N2 = n2;
      int N3 = n3;
      int N4 = n4; 
      return NT2_F77NAME(ilaenv)(&i, fname.c_str(), opts.c_str(), &N1, &N2, &N3, &N4,
                                 fname.size(), opts.size());
    }
  }
}
#endif
