/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_LAPACK_SV_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_LAPACK_SV_HPP_INCLUDED
#include <nt2/toolbox/linalg/details/utility/padding.hpp>
#include <nt2/toolbox/linalg/details/lapack/workspace.hpp>
#include <nt2/toolbox/linalg/details/lapack/sv.hpp>
#include <nt2/table.hpp>
#include <nt2/core/settings/forward/shape.hpp>
#include <nt2/core/settings/forward/storage_scheme.hpp>
#include <nt2/core/container/category.hpp>
#include <nt2/include/functions/extent.hpp>
#include <nt2/include/functions/size.hpp>
#include <nt2/include/functions/issquare.hpp>
#include <boost/preprocessor/cat.hpp>
#include <nt2/sdk/error/assert.hpp>

namespace nt2
{
  namespace ext
  {
    //catch all
    NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::sv_, tag::cpu_, 
                                (A)(SIZEA)(SHA)(STA)
                                (B)(SIZEB)(SHB)(STB), 
                                ((table_< unspecified_<A>, nt2::settings(SIZEA,SHA,STA) > ))
                                ((table_< unspecified_<B>, nt2::settings(SIZEB,SHB,STB) > ))
                                )
    {
      typedef long int result_type;
      typedef typename B::value_type value_type; 
      
      BOOST_FORCEINLINE result_type operator()(char const&, A&, B&)
      {
        BOOST_ASSERT_MSG(false, "The current lapack matrices forms/storages\
                                   is not supported by any sv lapack call");
        return 1; 
      }
    };
    
    NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::sv_, tag::cpu_, 
                                (A)(SIZEA)(SHA)(STA)
                                (B)(SIZEB)(SHB)(STB)
                                (C)(SIZEC)(SHC)(STC), 
                                ((table_< unspecified_<A>, nt2::settings(SIZEA,SHA,STA) > ))
                                ((table_< unspecified_<B>, nt2::settings(SIZEB,SHB,STB) > ))
                                ((table_< unspecified_<C>, nt2::settings(SIZEC,SHC,STC) > ))
                                )
    {
      typedef long int result_type;
      typedef typename B::value_type value_type; 
      
      BOOST_FORCEINLINE result_type operator()(char const&, A&, B&, C&)
      {
        BOOST_ASSERT_MSG(false, "The current lapack matrices forms/storages\
                                   is not supported by any sv lapack call");
      }
    };
    
    NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::sv_, tag::cpu_, 
                                (A)(SIZEA)
                                (B)(SIZEB), 
                                ((table_< floating_<A>, nt2::settings(SIZEA,nt2::rectangular_,nt2::conventional_) > ))
                                ((table_< floating_<B>, nt2::settings(SIZEB,nt2::rectangular_,nt2::conventional_) > ))
                                )
    {
      
      typedef long int result_type;
      typedef typename B::value_type value_type; 
      BOOST_FORCEINLINE result_type operator()(A & a, B & b)
      {
        typedef typename nt2::meta::make_container<nt2::tag::table_, long int, of_size_<-1>  >::type IPIV;
        IPIV ipiv(of_size(size(a, 1), 1));
        BOOST_ASSERT_MSG(is_square(a), "matrix a must be square in calling lapack sv routines"); 
        long int info;
        const long int n     = nt2::extent(a)[0];
        const long int nrhs  = nt2::extent(b)[1];
        const long int lda   = nt2::details::padding(a);
        const long int ldb   = nt2::details::padding(b);
        nt2::details::gesv (&n, &nrhs, a.begin(), &lda, ipiv.begin(), b.begin(), &ldb, &info);
        //          mcheck::LapackTest(__FILE__, __LINE__, "gesv", A, info); 
        return (info == 0);
        //        return nt2::sv(a, b, ipiv); 
      }
    };
    
    NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::sv_, tag::cpu_, 
                                (A)(SIZEA)
                                (B)(SIZEB)
                                (IPIV)(SIZEIPIV), 
                                ((table_< floating_<A>, nt2::settings(SIZEA,nt2::rectangular_,nt2::conventional_) > ))
                                ((table_< floating_<B>, nt2::settings(SIZEB,nt2::rectangular_,nt2::conventional_) > ))
                                ((table_< int32_<IPIV>, nt2::settings(SIZEIPIV,nt2::rectangular_,nt2::conventional_) > ))
                                )
    {
      typedef long int result_type;
      typedef typename B::value_type value_type; 
      BOOST_FORCEINLINE result_type operator()(A & a, B & b, IPIV& ipiv)
      {
        BOOST_ASSERT_MSG(is_square(a), "matrix a must be square in calling lapack sv routines"); 
        long int info;
        const long int n     = nt2::extent(a)[0];
        const long int nrhs  = nt2::extent(b)[1];
        const long int lda   = nt2::details::padding(a);
        const long int ldb   = nt2::details::padding(b);
        nt2::details::gesv (&n, &nrhs, a.begin(), &lda, ipiv.begin(), b.begin(), &ldb, &info);
        //          mcheck::LapackTest(__FILE__, __LINE__, "gesv", A, info); 
        return (info == 0);
      }
    };
  }
}
#endif
