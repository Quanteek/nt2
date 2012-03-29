/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_TOOLBOX_OPTIMIZE_FUNCTIONS_THRESHOLD_HPP_INCLUDED
#define NT2_TOOLBOX_OPTIMIZE_FUNCTIONS_THRESHOLD_HPP_INCLUDED
#include <nt2/sdk/complex/meta/as_real.hpp>
#include <nt2/include/constants/sqrteps.hpp>

namespace nt2
{

  template<typename T> class threshold
  {
  public :
    typedef T                                      float_t;
    typedef typename meta::as_real<float_t>::type   base_t;
    
    threshold(  base_t abs_ = Sqrteps<base_t>(),
                base_t rel_ = Sqrteps<base_t>(),
                base_t res_ = Sqrteps<base_t>(),
                size_t iter_ = 100)
      : absoluteTolerance(abs_), relativeTolerance(rel_), residualTolerance(res_), nbIterMax(iter_) {};
    
    base_t getAbsoluteTolerance()  const { return absoluteTolerance; }
    base_t getRelativeTolerance()  const { return relativeTolerance; }
    base_t getResidualTolerance()  const { return residualTolerance; }
    
    size_t  getNbIterMax()          const { return nbIterMax; }
    
    void setAbsoluteTolerance( base_t val )  { absoluteTolerance = val; }
    void setRelativeTolerance( base_t val )  { relativeTolerance = val; }
    void setResidualTolerance( base_t val )  { residualTolerance = val; }
    void setNbIterMax( size_t val )          { nbIterMax         = val; }
    
    ~threshold(){}; 
    
  protected :
    
    base_t absoluteTolerance;
    base_t relativeTolerance;
    base_t residualTolerance; 
    size_t  nbIterMax;
  };
  
}

#endif
