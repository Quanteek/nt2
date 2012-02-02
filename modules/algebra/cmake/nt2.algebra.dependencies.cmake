################################################################################
#         Copyright 2003 & onward LASMEA UMR 6602 CNRS/Univ. Clermont II
#         Copyright 2009 & onward LRI    UMR 8623 CNRS/Univ Paris Sud XI
#
#          Distributed under the Boost Software License, Version 1.0.
#                 See accompanying file LICENSE.txt or copy at
#                     http://www.boost.org/LICENSE_1_0.txt
################################################################################

find_package(BLAS)

set(BLAS_FOUND TRUE)
set(NT2_ALGEBRA_DEPENDENCIES_FOUND ${BLAS_FOUND}) 
set(NT2_ALGEBRA_DEPENDENCIES_LIBRARIES ${BLAS_LIBRARIES})

set(NT2_ALGEBRA_DEPENDENCIES_EXTRA core)