//==============================================================================
//         Copyright 2003 - 2012 LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#ifndef NT2_SDK_DETAILS_TYPE_ID_HPP_INCLUDED
#define NT2_SDK_DETAILS_TYPE_ID_HPP_INCLUDED

/**
* @file
* @brief Defines types to string conversion utility functions
**/

#if (__GNUC__ && __cplusplus && __GNUC__ >= 3)
//==============================================================================
// Includes abi::__cxa_demangle
//==============================================================================
#include <cxxabi.h>
#endif

#include <string>
#include <cstdlib>
#include <typeinfo>
#include <iostream>
#include <boost/type_traits/is_const.hpp>
#include <boost/type_traits/is_reference.hpp>
#include <boost/type_traits/remove_reference.hpp>

namespace nt2 {  namespace details
{
  // INTERNAL ONLY
  // demangle a type name retrieved through typeid()
  inline std::string demangle(const char* name)
  {
    #if(__GNUC__ && __cplusplus && __GNUC__ >= 3)
    std::size_t len;
    int         stat;

    char* realname = ::abi::__cxa_demangle(name,NULL,&len,&stat);

    if(realname != NULL)
    {
      std::string out(realname);
      ::free(realname);
      return out;
    }
    else
    {
      return std::string("?");
    }
    #else
    std::string out(name);
    return out;
    #endif
  }

  // INTERNAL ONLY
  // Stream some indentation on a output stream
  inline std::ostream& indent(std::ostream& os, size_t depth)
  {
    for(size_t i=0; i<depth; ++i)
      os << "    ";

    return os;
  }
} }

namespace nt2
{
  /**
  * @brief Type name demangling function
  *
  * For any given type @c T, returns a human readable string containing the fully
  * qualified name of @c T.
  *
  * @tparam T   Type to turn into a string
  * @return @c std::string containing the type of @c T
  *
  * @usage
  * @include type_id.cpp
  *
  * This examples output:
  *
  * @code
  * char [21]
  * float
  * std::vector<long*, std::allocator<long*> >
  * @endcode
  **/
  template<typename T> inline std::string type_id()
  {
    std::string s = details::demangle(typeid(T).name());
    if(boost::is_const<typename boost::remove_reference<T>::type>::value)
      s += " const";
    if(boost::is_reference<T>::value)
      s += "&";
    return s;
  }

  /**
  * @brief Type name demangling function
  *
  * For any given value @c x of type @c T, returns a human readable string
  * containing the fully qualified name of @c T.
  *
  * @param  x   Value to analyze
  * @return @c  std::string containing the type of @c x
  *
  * @usage
  * @include type_id.cpp
  *
  * This examples output:
  *
  * @code
  * char [21]
  * float
  * std::vector<long*, std::allocator<long*> >
  * @endcode
  **/
  template<typename T> inline std::string
  type_id ( const T&
#if defined(DOXYGEN_ONLY)
            x
#endif
          )
  {
    return type_id<T>();
  }

  /**
  *
  * @brief Type name streaming function
  *
  * For any given type @c T, displays a human readable string containing the fully
  * qualified name of @c T on the standard output. Formatting is applied on this
  * output so template types and other complex structures are properly displayed.
  *
  * @tparam T   Type to display
  *
  * @usage
  * @include display_type.cpp
  *
  * This examples outpus:
  *
  * @code
  * char [21]
  * float
  * std::vector<
  *              long*
  *              ,std::allocator<
  *                              long*
  *                              >
  *            >
  * @endcode
  **/
  template<typename T> inline void display_type()
  {
    std::string s = type_id<T>();

    size_t depth = 0;
    bool prevspace = true;
    for(std::string::const_iterator it = s.begin(); it != s.end(); ++it)
    {
      switch(*it)
      {
        case '<':
          depth++;
          std::cout << *it;
          std::cout << '\n';
          details::indent(std::cout, depth);
          prevspace = true;
          break;

        case '>':
          depth--;
          std::cout << '\n';
          details::indent(std::cout, depth);
          std::cout << *it;
          prevspace = false;
          break;

        case ',':
          std::cout << '\n';
          details::indent(std::cout, depth);
          std::cout << *it;
          prevspace = true;
          break;

        case ' ':
          if(!prevspace)
            std::cout << *it;
          break;

        default:
          std::cout << *it;
          prevspace = false;
      }
    }
    std::cout << std::endl;
  }

  /**
  * @brief Type name streaming function
  *
  * For any given value @c x of type @c T, displays a human readable string
  * containing the fully qualified name of @c T on the standard output.
  * Formatting is applied on this output so template types and other complex
  * structures are properly displayed.
  *
  * @param x   Value to display type from
  *
  * @usage
  *
  * @include display_type.cpp
  *
  * This examples outpus:
  *
  * @code
  * char [21]
  * float
  * std::vector<
  *              long*
  *              ,std::allocator<
  *                              long*
  *                              >
  *            >
  * @endcode
  **/
  template<typename T> inline void
  display_type( const T&
#if defined(DOXYGEN_ONLY)
                x
#endif
              )
  {
    return display_type<T>();
  }
}

#endif
