################################################################################
#         Copyright 2003 & onward LASMEA UMR 6602 CNRS/Univ. Clermont II
#         Copyright 2009 & onward LRI    UMR 8623 CNRS/Univ Paris Sud XI
#
#          Distributed under the Boost Software License, Version 1.0.
#                 See accompanying file LICENSE.txt or copy at
#                     http://www.boost.org/LICENSE_1_0.txt
################################################################################
################################################################################
# Find the Blas
################################################################################

find_package(BLAS QUIET)

if(BLAS_FOUND)
  set(NT2_BLAS_FOUND TRUE) 
  set(NT2_BLAS_LIBRARIES ${BLAS_LIBRARIES})
  set(NT2_BLAS_LINK_FLAGS ${BLAS_LINKER_FLAGS})
elseif(NOT BLAS_FOUND)
  
  # Intel MKL
  if(NT2_BLAS_VENDOR STREQUAL "Intel")
    if(NT2_ARCH_X86_64)
      find_library(NT2_MKL_LP64 NAMES mkl_intel_lp64 PATHS ${NT2_BLAS_ROOT})
      set(NT2_BLAS_LIBRARIES ${NT2_BLAS_LIBRARIES} ${NT2_MKL_LP64}) 
    elseif(NT2_ARCH_X86)
      find_library(NT2_MKL_32 NAMES mkl_intel PATHS ${NT2_BLAS_ROOT})
      set(NT2_BLAS_LIBRARIES ${NT2_BLAS_LIBRARIES} ${NT2_MKL_32}) 
    endif()


    #if(TRUE)
      nt2_module_tool(is_multicore RESULT_VARIABLE RESULT_VAR OUTPUT_QUIET)
      if(RESULT_VAR EQUAL -1)
        set(NT2_ARCH_MULTICORE FALSE)
      elseif(RESULT_VAR EQUAL 1)
        set(NT2_ARCH_MULTICORE FALSE)
      else()
        set(NT2_ARCH_MULTICORE TRUE)
      endif()
      message(${NT2_ARCH_MULTICORE}) 
    #else()
    #endif()
      find_library(NT2_MKL_CORE NAMES mkl_core PATHS ${NT2_BLAS_ROOT})
  else()
  endif()
  
  if(NOT NT2_BLAS_FOUND)
    message("-- Warning: If Blas is available on the target,")
    message("-- please specify the location and the vendor")
    message("-- by defining NT2_BLAS_ROOT and NT2_BLAS_VENDOR.")
    message("-- (\"Intel\")")
  endif()

endif()