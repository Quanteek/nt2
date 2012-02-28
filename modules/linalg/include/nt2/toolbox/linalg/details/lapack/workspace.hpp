/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_DETAILS_LAPACK_WORKSPACE_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_DETAILS_LAPACK_WORKSPACE_HPP_INCLUDED
                                                               
////////////////////////////////////////////////////////////////////////////////
/*! \file workspace.hpp
  \brief Use of LAPACK workspace.
  
  LAPACK routines use from one to 4  workspace arrays
  - work    is an array of the primary element type (d, s, z, c) the size of
            which can be determined by a call with parameter ldwork set to -1
  - rwork   is an array of the base primary element, (d, s, d, s)
  - iwork   is an array of long int
  - bwork   is an array of long int (logical)
  
  This class allows to reuse a workspace allocated once for all.
  
  Each interface (_itf.hh files) provides two kinds of C++ eqivalent calls to
  LAPACK routines :
  
  The first contains no workspace parameter and do all allocations/destructions
  internally,  and of course at each call.
  
  The second has a last parameter w which is a workspace reference that can be
  passed allocated or not. The routine first reallocate the workspace (which
  means to do nothing if the sizes already allocated are greater  than needs)
  While in scope the workspace grows and is less and less really resized.
  
  To optimize reallocation one can in a learning stage  declare a void
  workspace for each  used type of data used, to use it constantly and to
  issue calls to showstate at the end of the program. Then to change the
  workspace initial construction with the obtained sizes.
  
  The constructors of the factorization class of lpp::algebra also have a
  reference workspace optionnal parameter. These classes can use it instead
  of their own internal workspace. In fact,  all calls can share the same
  unique workspace.
  
  The interface has been written to be used with nt2 library, but it is
  independant from nt2 (http://nt2.sourceforge.net)
  std::vector has been choosen as storage,  but the scheme can be  easily
  adapted to put storage in any container class as soon as :
  - storage is contiguous, 
  - actual size and allocation are separate notions (i.e. reallocation is not
  done when actual size is bigger than needed).
  
*/

#include <vector>
#include <iostream>
#include <nt2/include/functions/real.hpp>
#include <nt2/include/functions/max.hpp>


namespace nt2
{
  namespace details
  {
    //////////////////////////////////////////////////////////////////////////
    //       This function is to cast result returned by
    //       LAPACK with the wrong type
    //      
    //////////////////////////////////////////////////////////////////////////

    template < class T > long int getval(const T & v)
    {
      // FORTRAN give here a float, double or complex of these,
      // but there is in fact a long int stored here.
      return static_cast<long int>(real(v)); 
      
    }
    /**
     *  This class allows to use easily LAPACK workspace features
     */
    
    template < class T > class workspace
    {
    public :
      typedef T                                 type_t;
      typedef typename meta::as_real<T>::type   base_t;
      // TO DO
      // Is vector the good container ? Have we better in our pockets ?
      // Must it be a policy of workspace ?
      typedef std::vector < type_t >               m_t;
      typedef std::vector < base_t >              mb_t;
      typedef std::vector < long int >            mi_t;
      /**
       * Sized constructor.
       * It allows to provide workspace memory from the beginning of the process
       */
      workspace(size_t wsiz, size_t rwsiz, size_t iwsiz, size_t bwsiz = 0):
        w(nt2::max(wsiz, size_t(1))), 
        rw(nt2::max(rwsiz, size_t(1))), 
        iw(nt2::max(iwsiz, size_t(1))),
        bw(nt2::max(bwsiz, size_t(1))),
        mquery(-1)
      {}
      /**
       * Trivial constructor.
       * Only the main work area get a necessary small provision.
       * The others will be allocated on demand.
       */
      workspace():w(1), rw(1), iw(1), bw(1), mquery(-1){ }
      ~workspace(){}
      //! to resize the main work area
      void resizew (size_t size){ w.reserve(size); } 
      //! to resize the real work area
      void resizerw(size_t size){rw.reserve(size); }    
      //! to resize the integer work area 
      void resizeiw(size_t size){iw.reserve(size); }  
      //! to resize the logical work area 
      void resizebw(size_t size){bw.reserve(size); }  
      
      bool isok(size_t wsiz, size_t rwsiz, size_t iwsiz, size_t bwsiz)const{ //! is the workspace sufficiently large
        return (wsiz >= w.size()) && (rwsiz >= rw.size())
          && (iwsiz >= iw.size()) && (bwsiz >= bw.size()); 
      }
      void realloc(size_t wsiz, size_t rwsiz, size_t iwsiz, size_t bwsiz = 0){
        if(wsiz >   w.size())  w.reserve(wsiz);
        if(rwsiz > rw.size()) rw.reserve(rwsiz);
        if(iwsiz > iw.size()) iw.reserve(iwsiz);
        if(bwsiz > bw.size()) bw.reserve(bwsiz);
      }
      type_t   * getw(){return &w[0]; }
      base_t   * getrw(){return &rw[0]; }
      long int * getiw(){return &iw[0]; }
      long int * getbw(){return &bw[0]; }
      const long int & neededsize(long int siz = 0){
        return (mldw =  siz ? siz : getval(w[0]));
      }
      const long int & neededisize(long int siz = 0){
        return mldiw = siz ? siz : iw[0];
      }
      const long int & neededrsize(long int siz = 0){
        return mldrw = siz ? siz : (long int)rw[0];
      }
      const long int & neededbsize(long int siz = 0){
        return mldbw = siz ? siz : bw[0];
      }
      const long int * query()const{return &mquery; }
      
      void showstate(char unit =  'B')const{
        size_t total = w.size()*sizeof(T)+rw.size()*sizeof(base_t)+iw.size()*sizeof(long int);
        size_t ext = w.capacity()*sizeof(T)+rw.capacity()*sizeof(base_t)+
          iw.capacity()*sizeof(long int)+bw.capacity()*sizeof(long int);
        std::cout << "workspace is using a total amount of: ";
        std::cout << total/factor(unit) << " " << name(unit) << std::endl;
        std::cout <<  w.size() << " in w,  " << std::endl;
        std::cout << rw.size() << " in rw, " << std::endl;
        std::cout << iw.size() << " in iw, " << std::endl;
        std::cout << bw.size() << " in bw, " << std::endl;
        std::cout << "Total capacity workspace up to now is  : ";
        std::cout << ext/factor(unit) << " " << name(unit) << std::endl; 
        std::cout <<  w.capacity() << " in w,  " << std::endl;
        std::cout << rw.capacity() << " in rw, " << std::endl;
        std::cout << iw.capacity() << " in iw, " << std::endl;
        std::cout << bw.capacity() << " in bw, " << std::endl;
      }
    private :
      m_t   w;
      mb_t rw;
      mi_t iw;
      mi_t bw;
      long int mldw, mldiw, mldrw, mldbw;
      const long int mquery;
      static inline size_t factor(char unit)
      {
        switch (unit){
        case 'M' :
        case 'm' :
          return(1024*1024); 
          break; 
        case 'B' :
        case 'b' :
        case 'O' :
        case 'o' :
          return 1;
          break; 
        case 'K' :
        case 'k' :
        default :
          return 1024; 
        }
        return 1024; 
      }
      static inline std::string name(char unit)
      {
        switch (unit){
        case 'M' :
        case 'm' :
          return("Megabytes"); 
          break; 
        case 'B' :
        case 'b' :
        case 'O' :
        case 'o' :
          return("Bytes"); 
          break; 
        case 'K' :
        case 'k' :
        default :
          return("Kilobytes"); 
        }
        return("Kilobytes"); 
      }
    };
  }    
}
#endif
