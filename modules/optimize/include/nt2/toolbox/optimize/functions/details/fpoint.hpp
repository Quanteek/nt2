/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_OPTIMIZE_FUNCTIONS_DETAILS_FPOINT_HPP_INCLUDED
#define NT2_TOOLBOX_OPTIMIZE_FUNCTIONS_DETAILS_FPOINT_HPP_INCLUDED
#include <nt2/include/functions/optimizable.hpp>
#include <nt2/include/constants/one.hpp>
#include <nt2/include/functions/abs.hpp>
#include <nt2/include/functions/first_index.hpp>
#include <nt2/include/functions/last_index.hpp>   
#include <nt2/include/functions/numel.hpp>
#include <nt2/table.hpp>

namespace nt2
{
  template < class T, class FLOAT, class FUNC> class fpoint
  {
    typedef T                                  matrix_t;
    typedef FLOAT                               float_t;
    typedef nt2::table<float_t>                 table_t;
    typedef FUNC                                 func_t;
    typedef eval<matrix_t,float_t,func_t>       eval_t;
    typedef fpoint<matrix_t,float_t,func_t>      self_t;
    
    matrix_t& x;                                // A point in the function's domain.
    float_t fval;                               // Function value at the point.
    const func_t & fproc;                       // Procedure to compute that value
    const bool free_x_on_destructing;           // The flag telling if this fpoint 
                                                // "owns" x, and has to dispose of
                                                // its dynamic memory on destruction.
      
      
  public:
    fpoint(matrix_t& b, const func_t & f)
      : x(b),
        fproc(f),
        free_x_on_destructing(false)
    {
      eval_t::DoEval(fval, f, b);
    }
    
    ~fpoint(){
      if( free_x_on_destructing ) delete &x;
    };
    
    float_t f() const { return fval; }
    
    fpoint(const fpoint& fp) :
      x(*(new matrix_t(fp.x))), 
      fval(fp.fval),
      fproc(fp.fproc),
      free_x_on_destructing(true){}
    
    fpoint& operator = (const fpoint& fp){{
        if (&fp != this){
          x = fp.x;
          fval = fp.fval;
        }
        return *this;
      }
    };
    
    float_t fiddle_around(const matrix_t& h);
    // Examine the function in the
    // neighborhood of the current point.
    // h defines the radius of the region
    
    // Proceed in the direction the function
    // seems to decline
    template < class FPOINT > 
    friend void update_in_direction(FPOINT & from, FPOINT & to);
    
    // Decide whether the region embracing
    // the local min is small enough
    bool isstep_relatively_small(const matrix_t& h, const float_t tau)
    {
      for(ptrdiff_t i=first_index<0>(x); i <= last_index<0>(x) ; i++)
        {
          if (float_t(h(i)) >= float_t(tau*(One<float_t>()+nt2::abs(x(i))))) return false; 
        }
      return true;
      //TODO nt2::all(h < tau*(One<float_t>()+nt2::abs(x))); 
    }
  };
  
  /*
   * Examine the function f in the vicinity of the current point x
   * by making tentative steps fro/back along each coordinate.
   * Should the function decrease, x is updated to pinpoint thus
   * found new local min.
   * The procedure returns the minimal function value found in
   * the region.
   *
   */
  template < class T, class FLOAT, class FUNC  > 
  FLOAT fpoint<T, FLOAT, FUNC >::fiddle_around(const matrix_t& h)
  {
    // Perform a step along each coordinate
    for(size_t i = nt2::first_index<0>(x); i <= nt2::last_index<0>(x); ++i)
      { 
        const float_t hi = h(i); 
        const float_t xi_old = x(i);          // Old value of x[i]
        float_t fnew;
        x(i) = xi_old + hi;
        eval_t::DoEval(fnew, fproc, x); 
        if ( fnew < fval )// Step caused f to decrease, OK
          {
            fval = fnew; return fval;  
          }
        x(i) = xi_old - hi;
        eval_t::DoEval(fnew, fproc, x); 
        if (fnew < fval )// Step caused f to decrease, OK
          {
            fval = fnew; return fval; 
          }
        x(i) = xi_old; // No function decline has been found along this coord, back up
      }
    return fval;
  }                                                
  
  // Proceed in the direction the function
  // seems to decline
  // to_new = (to - from) + to
  // from = to (before modification)
  template < class FPOINT  > 
  void update_in_direction(FPOINT& from, FPOINT& to)
  {
    typedef typename FPOINT::float_t float_t; 
    for(size_t i = nt2::first_index<0>(from.x); i <= nt2::last_index<0>(from.x); ++i)
      {
        const float_t t = to.x(i);
        to.x(i)  += (t - from.x(i));
        from.x(i) = t;
      }
    from.fval = to.fval;
    FPOINT::eval_t::DoEval(to.fval, to.fproc, to.x); 
  }
}

#endif
