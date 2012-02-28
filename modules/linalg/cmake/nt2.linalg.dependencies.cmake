################################################################################
#         Copyright 2003 & onward LASMEA UMR 6602 CNRS/Univ. Clermont II
#         Copyright 2009 & onward LRI    UMR 8623 CNRS/Univ Paris Sud XI
#
#          Distributed under the Boost Software License, Version 1.0.
#                 See accompanying file LICENSE.txt or copy at
#                     http://www.boost.org/LICENSE_1_0.txt
################################################################################

include(nt2.blas)

set(NT2_LINALG_DEPENDENCIES_FOUND ${NT2_BLAS_FOUND}) 
set(NT2_LINALG_DEPENDENCIES_LIBRARIES ${NT2_BLAS_LIBRARIES})
set(NT2_LINALG_LINK_FLAGS ${NT2_BLAS_LINK_FLAGS})

set(NT2_LINALG_DEPENDENCIES_EXTRA sdk core openmp)