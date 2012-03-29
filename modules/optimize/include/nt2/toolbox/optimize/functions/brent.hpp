/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_OPTIMIZE_FUNCTIONS_BRENT_HPP_INCLUDED
#define NT2_TOOLBOX_OPTIMIZE_FUNCTIONS_BRENT_HPP_INCLUDED
#include <nt2/include/functions/threshold.hpp>
#include <nt2/include/functions/optimizable.hpp>
#include <nt2/include/functions/average.hpp>
#include <nt2/include/functions/sign.hpp>
#include <nt2/include/functions/abs.hpp>
#include <nt2/include/constants/half.hpp>
#include <nt2/include/constants/two.hpp>
#include <nt2/include/constants/zero.hpp>
#include <nt2/include/constants/cgold.hpp>
#include <nt2/include/constants/nan.hpp>


namespace nt2
{
  template<typename FLOAT> class brent : public threshold<FLOAT> 
  {
  public :
    typedef FLOAT float_t;
    
    brent(size_t itmax_ = 100):itmax(itmax_), iter(0){} 
    ~brent() {}
    
    template < class FUNC > void optimize( const FUNC & f,
                                           float_t low_,
                                           float_t init_,
                                           float_t up_);
    
    size_t  getNbIteration()      const { return iter;   }
    float_t getMinimum()          const { return valmin; }
    float_t getMinimumPosition()  const { return xmin;   } 
    float_t getLowerBound()       const { return ax;     }
    float_t getUpperBound()       const { return cx;     }
    float_t getInitialPoint()     const { return bx;     }
    
  private:
    
    void shift(float_t& a, float_t& b,
               float_t& c, float_t d)
    {
      a = b; b = c; c = d;
    }
    
    const size_t itmax;
    size_t        iter; 
    float_t ax, bx, cx;
    float_t       xmin;
    float_t     valmin;
  };    
  
  
  template<typename FLOAT>
  template<typename FUNC> 
  inline void brent<FLOAT>::optimize( const FUNC& f, float_t low_, float_t init_, float_t up_ )
  {
    ax = low_; bx = init_;  cx = up_; iter = 0; 
    float_t d = Nan<float_t>() ,etemp,fu,p1,q,r,u;
    float_t e = Zero<float_t>();
    float_t a = (ax < cx ? ax : cx);
    float_t b = (ax > cx ? ax : cx);
    float_t v,w,x; 
    x=w=v=bx;
    float_t fx, fw, fv; 
    eval<float_t,float_t,FUNC>::DoEval(fx,f,x); 
    fw=fv=fx;
    for(iter=1;iter<=itmax;++iter) 
      {
        float_t tol = this->getAbsoluteTolerance(); //outside the loop ?
        float_t xm=nt2::average(a, b);
        float_t tol1=tol*nt2::abs(x)+Sqrteps<float_t>();
        float_t tol2=Two<float_t>()*tol1;
        if (nt2::abs(x-xm) <= (tol2-Half<float_t>()*(b-a))) 
          {
            xmin   = x;
            valmin = fx;
            return;
          }
        if (nt2::abs(e) > tol1) 
          {
            r=(x-w)*(fx-fv);
            q=(x-v)*(fx-fw);
            p1=(x-v)*q-(x-w)*r;
            q=Two<float_t>()*(q-r);
            if (q > Zero<float_t>()) p1 = -p1;
            q=nt2::abs(q);
            etemp=e;
            e=d;
            
            if (nt2::abs(p1) >= nt2::abs(Half<float_t>()*q*etemp) ||
                (p1 <= q*(a-x)) ||
                (p1 >= q*(b-x))) 
              {
                d=Cgold<float_t>()*(e=(x >= xm ? a-x : b-x));
              }
            else 
              {
                d=p1/q;
                u=x+d;
                if (u-a < tol2 || b-u < tol2) d=tol1 * nt2::sign(xm-x);
              }
          } 
        else 
          {
            d=Cgold<float_t>()*(e=((x >= xm) ? a-x : b-x));
          }
        
        u=(nt2::abs(d) >= tol1 ? x+d : x + tol1 *  nt2::sign(d));
        //                    fu= f(u);
        eval< float_t,float_t,FUNC>::DoEval(fu,f,u); 
        //      std::cout << "u " << u <<  "  fu " <<  fu << std::endl; 
        if (fu <= fx) 
          {
            if (u >= x) a=x; else b=x;
            shift(v,w,x,u);
            shift(fv,fw,fx,fu);
          } 
        else 
          {
            if (u < x) a=u; else b=u;
            if (fu <= fw || w == x)
              {
                v=w;
                w=u;
                fv=fw;
                fw=fu;
              } 
            else if (fu <= fv || v == x || v == w) 
              {
                v=u;
                fv=fu;
              }
          }
      }
    BOOST_ASSERT_MSG(false, "Brent has not converged"); 
    xmin   = x;
    valmin = fx;
  }
}

#endif
