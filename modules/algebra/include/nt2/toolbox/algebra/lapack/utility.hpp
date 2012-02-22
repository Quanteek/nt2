/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_ALGEBRA_LAPACK_UTILITY_HPP_INCLUDED
#define NT2_TOOLBOX_ALGEBRA_LAPACK_UTILITY_HPP_INCLUDED
#include <string>
#include <nt2/toolbox/algebra/blas/f77_wrapper.hpp > 

namespace nt2
{
  namespace lpp
  {
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
  
#include "lapackworkspace.hh"
  
  ////////////////////////////////////////////////////////////////
  // csrot  rotation plane à deux vecteurs complexes ??
  // zdrot  rotation plane à deux vecteurs complexes ??
  // csrscl multiplication vecteur 1/a,  inutile
  // zdrscl multiplication vecteur 1/a,  inutile
  // xerbla gestion d'erreurs
  // second temps en seconde,  inutile
  // ladiv  à faire
  ////////////////////////////////////////////////////////////////
  // interfacer lamch avec les constantes
  ////////////////////////////////////////////////////////////////
  
  namespace lpp
  {
    inline char lower(const char c){
      return (c >= 'A' && c <= 'Z') ? c +('a'-'A') : c; 
    }
    inline char upper(const char c){
      return (c >= 'a' && c <= 'z') ? c -('a'-'A') : c; 
    }
    
  }
  // /////////////////////////////////////////////////////////////////////////////
  //  End of lpp namespace
  // /////////////////////////////////////////////////////////////////////////////
}
#endif
