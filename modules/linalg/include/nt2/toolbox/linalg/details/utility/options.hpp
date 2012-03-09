/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_LINALG_DETAILS_UTILITY_OPTIONS_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_DETAILS_UTILITY_OPTIONS_HPP_INCLUDED

namespace nt2
{
  namespace details
  {

//     template < class T > struct lapack {
//       typedef typename ttt::type<T>::base_t                                  base_t; 
//       static const bool is_lapack = ttt::type<T>::is_float || ttt::type<base_t>::is_float; 
//     };
    
//     const nc::matrix < char >  laoptions_abc = nt2::colon('A', 'Z');
    
//     static inline const char & opt(char c){
//       static char laoptions_abc[] = {65, 66, 67, 68, 69,
//                        70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
//                        80, 81, 82, 83, 84, 85, 86, 87, 88, 89,
//                        90, 91}; 
//       return laoptions_abc[c-65]; 
//     }
    
    class options{
    public:
      options(bool trans =  false):
        lt(false),
        ut(false),
        uhess(false),
        sym(false),
        posdef(false),
        rect(true),
        transa(trans),
        unidiag(false),
        general(true)
      {}
      
//      template < class XPR > 
//       options(const expression < XPR > &a, bool trans)//:
//        lt(istril(a)),
//         ut(!lt && istriu(a)),
//         uhess(!ut && isuhess(a)),
//         sym(!ut && issym(a)),
//         posdef(ut && ishermposdiag(a)),
//         rect(a.width()!= a.height()),
//         transa(trans)
//       {

//         lt = (istril(a));
//        ut = (!lt && istriu(a));
//         uhess = (!ut && isuhess(a));
//         sym = (!ut && issym(a));
//         posdef = (ut && ishermposdiag(a));
//         rect = (a.width()!= a.height());
//         transa = (trans); 
 

//       }
      
      void setLT(bool b){
        lt = b;
        if (lt){
          general = posdef = sym = uhess = ut = false;
        }
      }
      void setUT(bool b){
        ut = b;
        if (ut){
          general = posdef = sym = uhess = lt = false;
        }
      }
      void setUHESS(bool b){
        uhess = b;
        if (ut){
          general = rect = posdef = sym = lt = ut = false;
        }
      }
      void setSYM(bool b){
        sym = b;
        if (sym){
          general = uhess = lt = ut = false;
        }
      }
     void setRECT(bool b){
        rect = b;
        if (rect){
          sym = posdef = uhess = lt = ut = false;
          general =  true; 
        }
      }
      
      void setPOSDEF(bool b){
        posdef = b; 
        if (posdef){
          general = uhess = lt = ut = rect = false;
        }
      }
      
      void setUNIDIAG(bool b){
        unidiag = b;
      }
      
      void setTRANSA(bool b){
        transa = b;
      }
      bool getUT()const{return ut; }
      bool getLT()const{return lt; }
      bool getUHESS()const{return uhess; }
      bool getSYM()const{return sym; }
      bool getPOSDEF()const{return posdef; }
      bool getTRANSA()const{return transa; }
      bool getUNIDIAG()const{return unidiag; }
      bool getRECT()const{return rect; }
      
      const char& getUplo()const{return opt(ut?'U':'L'); }
      
      template < class T >
      const char& getTrans(const T & )const{
        return opt(transa ? (!is_real(T(1))? 'C':'T') :'N');
      }
      
      const char& getUni()const{return opt(unidiag?'U':'N'); }
       
      ~options(){}; 
    private:
      static inline const char & opt(char c){
        static char laoptions_abc[] = {65, 66, 67, 68, 69,
                                       70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
                                       80, 81, 82, 83, 84, 85, 86, 87, 88, 89,
                                       90, 91}; 
        return laoptions_abc[c-65]; 
    }
      
      bool lt;
      bool ut;
      bool uhess;
      bool sym;
      bool posdef;
      bool rect;
      bool transa;
      bool unidiag;
      bool general; 
    }; 
      
    
  }
}
#endif
