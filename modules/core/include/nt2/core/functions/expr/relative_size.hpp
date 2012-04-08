//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#ifndef NT2_CORE_FUNCTIONS_EXPR_RELATIVE_SIZE_HPP_INCLUDED
#define NT2_CORE_FUNCTIONS_EXPR_RELATIVE_SIZE_HPP_INCLUDED

#include <nt2/core/container/dsl.hpp>
#include <nt2/core/container/category.hpp>
#include <nt2/include/functions/numel.hpp>

namespace nt2 { namespace ext
{
  //============================================================================
  // Compute indexing size using any expression in the non 1D cases
  //============================================================================
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::relative_size_, tag::cpu_
                            , (Idx)(Tag)(Arity)(Size)(Current)(Dims)
                            , ((expr_ < unspecified_<Idx>
                                      , Tag
                                      , Arity
                                      >
                              ))
                              (fusion_sequence_<Size>)
                              (mpl_integral_< scalar_< integer_<Current> > >)
                              (mpl_integral_< scalar_< integer_<Dims> > >)
                            )
  {
    typedef typename meta::call<tag::numel_(Idx const&)>::type  result_type;
      
    BOOST_DISPATCH_FORCE_INLINE result_type
    operator()(const Idx& idx, const Size&, const Current&, const Dims& ) const
    {
      return nt2::numel( idx );
    }
  };

  //============================================================================
  // Compute indexing size using _
  //============================================================================
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::relative_size_, tag::cpu_
                            , (Idx)(Arity)(Size)(Current)(Dims)
                            , ((expr_ < colon_< Idx >
                                      , nt2::tag::terminal_
                                      , Arity
                                      >
                              ))
                              (fusion_sequence_<Size>)
                              (mpl_integral_< scalar_< integer_<Current> > >)
                              (mpl_integral_< scalar_< integer_<Dims> > >)
                            )  
  {
    //==========================================================================
    // If _ is the last indexer, return the slice of all remaining sizes
    //==========================================================================
    typedef meta::call<tag::numel_( Size, boost::mpl::int_<Current::value> )> true_type;

    //==========================================================================
    // Else, return current size
    //==========================================================================
    typedef boost::fusion::result_of::at_c<Size const, Current::value>  false_type;

    //==========================================================================
    // Is this _ the final one in the indexers list ?
    //==========================================================================
    typedef boost::mpl::bool_<(Current::value == Dims::value-1)>      is_final;

    //==========================================================================
    // Select result_type accordingly
    //==========================================================================
    typedef typename  boost::mpl::
                      eval_if<is_final, true_type, false_type>::type  result_type;
    
    //==========================================================================
    // Main function entry point
    //==========================================================================
    BOOST_DISPATCH_FORCE_INLINE result_type
    operator()(const Idx&, const Size& sz, const Current& c, const Dims&) const
    {
      return eval(sz,c,is_final());
    }

    //==========================================================================
    // Dispatched implementation
    //==========================================================================
    BOOST_DISPATCH_FORCE_INLINE result_type
    eval(const Size& sz, const Current&, boost::mpl::true_ const&) const
    {
      return numel(sz, boost::mpl::int_<Current::value>());
    }

    BOOST_DISPATCH_FORCE_INLINE result_type
    eval(const Size& sz, const Current&, boost::mpl::false_ const&) const
    {
      return boost::fusion::at_c<Current::value>(sz);
    }
  };
} }

#endif
