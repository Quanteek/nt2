/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_OPTIMIZE_FUNCTIONS_OPTIMIZABLE_HPP_INCLUDED
#define NT2_TOOLBOX_OPTIMIZE_FUNCTIONS_OPTIMIZABLE_HPP_INCLUDED

namespace nt2
{
  
  template<typename T_IN, typename T_OUT = T_IN>
  struct optimizableProc
  {
    typedef T_IN   tin_t;
    typedef T_OUT  tout_t; 
    static const size_t arity = 2;  
    virtual void operator()( tout_t& fx, const tin_t& x) const = 0;
    virtual ~optimizableProc(){}
  };
  
  template<typename T_IN, typename T_OUT = T_IN>
  struct optimizableFunc
  {
    typedef T_IN   tin_t;
    typedef T_OUT  tout_t; 
    static const size_t arity = 1;  
    virtual tout_t operator()( const tin_t& x) const = 0;
    virtual ~optimizableFunc(){}
  };
  
  template <  class T_IN, class T_OUT >
  struct funTypes{
    typedef T_OUT (*func_t)(const T_IN&);
    typedef void  (*proc_t)( T_OUT&, const T_IN&);
  };
  
  template <  class T_IN, class T_OUT, class T_AUX = T_IN>
  struct jfunTypes{
    typedef T_OUT (*jfunc_t)(const T_IN&, const T_AUX & );
    typedef void  (*jproc_t)( T_OUT&, T_IN&, const T_AUX & );
  };
  
  template < class T_IN, class T_OUT, class F, class T_AUX = T_IN>
  struct evalf {
    typedef typename funTypes < T_IN, T_OUT > ::proc_t proc_t; 
    static void DoEval( T_OUT & y, const proc_t & f, const T_IN &x)
    {
      f(y, x);
    }
    typedef typename jfunTypes < T_IN, T_OUT, T_AUX> ::jproc_t jproc_t; 
    static void DoEval( T_OUT & y, const jproc_t & f, T_IN &x, const T_AUX &aux)
    {
      f(y, x, aux);
    }
  }; 
  
  
  
  template < class T_IN, class T_OUT, class T_AUX >
  struct evalf < T_IN,  T_OUT, typename funTypes<T_IN,T_OUT>::func_t, T_AUX> {
    typedef typename funTypes < T_IN, T_OUT > ::func_t func_t; 
    static void DoEval( T_OUT & y, const func_t& f, const T_IN & x)
    {
      y = f(x); 
    }
    typedef typename jfunTypes < T_IN, T_OUT, T_AUX> ::jproc_t jproc_t; 
    static void DoEval( T_OUT & y, const jproc_t & f, const T_IN &x, const T_AUX &aux)
    {
      y = f(x, aux);
    }
  };
  
  
  template < class T_IN, class T_OUT, class F, class T_AUX = T_IN, size_t ARITY = F::arity>
  struct evalmf
  {
    static void DoEval( T_OUT & y, const F &f, const T_IN x)
    {
      y = f(x); 
    }
    static void DoEval( T_OUT & y, const F &f, const T_IN x, const T_AUX aux)
    {
      y = f(x, aux); 
    }
  }; 
  
  template < class T_IN, class T_OUT, class F,class T_AUX >
  struct evalmf< T_IN,  T_OUT, F, T_AUX, 2>
  {
    static void DoEval( T_OUT & y, F f, const T_IN x)
    {
      f(y, x); 
    }
    static void DoEval( T_OUT & y,  F f, T_IN x, const T_AUX aux)
    {
      f(y, x, aux); 
    }
  }; 
  
  template <  class T_IN, class T_OUT,  class F,  class T_AUX = T_IN, bool ISMF = boost::is_class<F>::value >
  struct eval
  {
    static void DoEval( T_OUT & y,  const F &f, const T_IN & x)
    {
      evalmf< T_IN,  T_OUT, F >::DoEval(y, f, x); 
    }
    static void DoEval( T_OUT & y,  const F &f, const T_IN & x, const T_AUX aux)
    {
      evalmf< T_IN,  T_OUT, F, T_AUX >::DoEval(y, f, x, aux); 
    }
    
  };
  
  template <  class T_IN, class T_OUT,  class F,  class T_AUX >
  struct eval < T_IN,  T_OUT, F, T_AUX, false > 
  {
    static void DoEval( T_OUT & y, const F& f, const T_IN & x)
    {
      evalf < T_IN, T_OUT, F > ::DoEval(y, f, x); 
    }
    static void DoEval( T_OUT & y,  const F &f, const T_IN & x, const T_AUX aux)
    {
      evalf< T_IN,  T_OUT, T_AUX, F >::DoEval(y, f, x, aux); 
    }
    
  };
  
}

#endif
