//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#ifndef NT2_SDK_MEMORY_DETAILS_ID_CONTAINER_HPP_INCLUDED
#define NT2_SDK_MEMORY_DETAILS_ID_CONTAINER_HPP_INCLUDED

#include <nt2/sdk/memory/details/container_base.hpp>

namespace nt2 { namespace memory
{
  //============================================================================
  // When given a non-null ID, containers share their block as a static member.
  // Every related methods and values are therefore static.
  //============================================================================
  template<class Tag, class ID, class T, class S>
  struct container : nt2::details::container_base<Tag,T,S>
  {
    typedef nt2::details::container_base<Tag,T,S>         parent;
    typedef typename parent::block_t                      block_t;
    typedef typename parent::allocator_type               allocator_type;
    typedef typename parent::value_type                   value_type;
    typedef typename parent::sizes_type                   sizes_type;
    typedef typename parent::extent_type                  extent_type;
    typedef typename parent::size_type                    size_type;
    typedef typename parent::difference_type              difference_type;
    typedef typename parent::reference                    reference;
    typedef typename parent::const_reference              const_reference;
    typedef typename parent::iterator                     iterator;
    typedef typename parent::const_iterator               const_iterator;
    typedef typename parent::pointer                      pointer;
    typedef typename parent::const_pointer                const_pointer;
    typedef typename parent::specific_data_type           specific_data_type;

    //==========================================================================
    // First constructor call is given priority over the others
    // Default constructor never throw nor assert as multiple instance can
    // coexist
    //==========================================================================
    container( allocator_type const& a = allocator_type())
    {
      if(!status_)
      {
        typename parent::block_t that(a);
        block_.swap(that);  // bleh swap :s
        parent::init( block_,sizes_
                    , typename parent::require_static_init()
                    );
        status_ = true;
      }
    }

    template<class Size>
    container ( Size const& sz
              , allocator_type const& a = allocator_type()
              )
    {
      if(!status_)
      {
        typename parent::block_t that(a);
        block_.swap(that);
        parent::init(block_,sz);
        sizes_  = sz;
        status_ = true;
      }
    }

    //==========================================================================
    // If I construct some container<ID> from another one, it is already in
    //==========================================================================
    container( container const& ) {}

    //==========================================================================
    // If I assign some container<ID> to another one, it is already in
    //==========================================================================
    container& operator=( container const& ) { return *this; }

    //==========================================================================
    // If I swap some container<ID> with another one, it is already in
    //==========================================================================
    void swap( container & ) {}

    //==========================================================================
    // Element access from instance
    //==========================================================================
    template<class Position> BOOST_FORCEINLINE
    reference operator[]( Position const& pos )
    {
      return parent::access(pos,block_,sizes_);
    }

    template<class Position> BOOST_FORCEINLINE
    const_reference operator[]( Position const& pos ) const
    {
      return parent::access(pos,block_,sizes_);
    }

    //==========================================================================
    // Return the ith index
    //==========================================================================
    BOOST_FORCEINLINE pointer
    get( difference_type i )
    {
      return block_.get(i);
    }

    BOOST_FORCEINLINE const_pointer
    get( difference_type i ) const
    {
      return block_.get(i);
    }

    //==========================================================================
    // Size of the container
    //==========================================================================
    static BOOST_FORCEINLINE extent_type extent() { return sizes_; }

    static BOOST_FORCEINLINE size_type size()       { return block_.size();       }

    static BOOST_FORCEINLINE bool empty() { return block_.empty(); }

    static BOOST_FORCEINLINE specific_data_type get_spec_data() { return specific_data_; }

    //==========================================================================
    /*!
     * Return the begin of the raw memory
     */
    //==========================================================================
    BOOST_FORCEINLINE pointer         raw()          { return block_.raw();    }
    BOOST_FORCEINLINE const_pointer   raw() const    { return block_.raw();    }
    BOOST_FORCEINLINE iterator        begin()        { return block_.begin();  }
    BOOST_FORCEINLINE const_iterator  begin() const  { return block_.begin();  }
    BOOST_FORCEINLINE iterator        end()          { return block_.end();    }
    BOOST_FORCEINLINE const_iterator  end() const    { return block_.end();    }

    //==========================================================================
    /*!
     * Return the number of physical element on the leading dimension
     */
    //==========================================================================
    static BOOST_FORCEINLINE size_type leading_size()
    {
      return parent::leading_size(sizes_);
    }

    //==========================================================================
    // Resize of the container
    //==========================================================================
    template<class Size> static BOOST_FORCEINLINE void resize( Size const& szs )
    {
      parent::resize(block_,szs,sizes_);
    }

    private:
    static bool                                 status_;
    static typename parent::block_t             block_;
    static typename parent::sizes_type          sizes_;
    static typename parent::specific_data_type  specific_data_;
  };

  //========================================================================
  // By default, block is non initialized
  //========================================================================
  template<class Tag, class ID, class T, class S>
  typename nt2::details::container_base<Tag,T,S>::block_t
  container<Tag, ID, T, S>::block_;

  //========================================================================
  // Set status to not constructed first
  //========================================================================
  template<class Tag, class ID, class T, class S>
  bool container<Tag, ID, T, S>::status_ = false;

  //========================================================================
  // Set size to default value
  //========================================================================
  template<class Tag, class ID, class T, class S>
  typename nt2::details::container_base<Tag,T,S>::sizes_type
  container<Tag, ID, T, S>::sizes_;

  //========================================================================
  // Set specific_data to default value
  //========================================================================
  template<class Tag, class ID, class T, class S>
  typename nt2::details::container_base<Tag,T,S>::specific_data_type
  container<Tag, ID, T, S>::specific_data_;

} }

#endif
