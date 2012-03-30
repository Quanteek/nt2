/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_OPTIMIZE_FUNCTIONS_HJMIN_HPP_INCLUDED
#define NT2_TOOLBOX_OPTIMIZE_FUNCTIONS_HJMIN_HPP_INCLUDED
#include <nt2/include/functions/threshold.hpp>
#include <nt2/include/functions/optimizable.hpp>
#include <nt2/toolbox/optimize/functions/details/fpoint.hpp>
#include <nt2/include/functions/threshold.hpp>
#include <nt2/include/functions/repnum.hpp>
#include <nt2/include/constants/oneoten.hpp>
#include <nt2/table.hpp>

namespace nt2
{
  template<typename T, typename FLOAT = typename T::value_type>
  struct hjmin : public threshold<FLOAT> 
  {
  public :
    typedef FLOAT                                   float_t;
    typedef T                                      matrix_t;
    typedef nt2::table<float_t>                     table_t;
    
    hjmin() : threshold<float_t>(), step_reduce_factor(Oneoten<float_t>()) {}
    
    ~hjmin() {}
    
    template < class FUNC >
    void optimize( const FUNC& crit,
                   matrix_t& aa,               //unknowns initialization
                   matrix_t& hh   
                   ); 
    template < class FUNC >
    void optimize( const FUNC& crit,
                   matrix_t& aa,               //unknowns initialization
                   const float_t h0   
                   ); 
    void setStepReduceFactor(float_t rf_){step_reduce_factor = rf_;      } 
    size_t               getNbIteration()      const { return iterdone;  }
    float_t              getMinimum()          const { return valmin;    }
    matrix_t             getMinimumPosition()  const { return *xmin;     } 
    
  private :
    float_t step_reduce_factor;
    size_t            iterdone; 
    matrix_t *            xmin;
    float_t             valmin; 
  }; 
    
  template<typename T, typename FLOAT>
  template < class FUNC > 
  void hjmin<T, FLOAT>::optimize( const FUNC& crit, T &aa, T &h)
  {
    xmin =  &aa;
    iterdone = 0; 
    const float_t tau = this->getRelativeTolerance();                 // Termination criterion.
    const float_t threshold = this->getAbsoluteTolerance();
    fpoint < matrix_t, float_t, FUNC> pmin(aa,crit);                  // Point of min.
    fpoint < matrix_t, float_t, FUNC> pbase(pmin);                    // Base point.
    for(size_t i = 1; i <= this->getNbIterMax() ; ++i)                // Main iteration loop.
      {                                                               // pmin is the approximation to min so far.
        iterdone = i;
        if( pbase.fiddle_around(h) < pmin.f() - threshold )
          {                                                           // Function value dropped significantly
            do                                                        // from pmin to the point pbase
              update_in_direction(pmin,pbase);                        // Keep going in the same direction
            while( pbase.fiddle_around(h) < pmin.f() - threshold );   // while it works
            pbase = pmin;                                             // Save the best approx found
          }
        else                                                          // Function didn't fall significantly
          {                                                           // upon wandering around pbas
            h = h*step_reduce_factor; 
            if( pbase.isstep_relatively_small(h,tau) )
              {
                valmin = pmin.f();
                return;
              }
          }
      }
    BOOST_ASSERT_MSG(false, "hjmin was not convergent"); 
  }
  
  // The same as above with the only difference
  // initial steps are given to be the same
  // along every direction. The final steps
  // aren't reported back though
  template<typename T, typename FLOAT> template < class FUNC > 
  void hjmin<T, FLOAT>::optimize( const FUNC& crit, T& aa, const FLOAT h0)
  {
    matrix_t h = nt2::repnum(h0, nt2::of_size(1, numel(aa)));
    return optimize(crit, aa, aa);
  }  
}

#endif
