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
include(nt2.info)

set(NT2_BLAS_FOUND FALSE)

  # Intel MKL
  if(NT2_BLAS_VENDOR STREQUAL "Intel")
  
    if(NT2_ARCH_X86_64)
      find_library(NT2_MKL_LP64 NAMES mkl_intel_lp64 mkl_intel_lp64_dll 
                   PATHS /opt/intel/mkl/lib/ ${NT2_BLAS_ROOT}
                   PATH_SUFFIXES ia32 intel64
                   )
      set(NT2_BLAS_LIBRARIES ${NT2_BLAS_LIBRARIES} ${NT2_MKL_LP64}) 
    elseif(NT2_ARCH_X86)
      find_library(NT2_MKL_32 NAMES mkl_intel mkl_intel_c_dll
                   PATHS /opt/intel/mkl/lib/ ${NT2_BLAS_ROOT}
                   PATH_SUFFIXES ia32 intel64
                   )
      set(NT2_BLAS_LIBRARIES ${NT2_BLAS_LIBRARIES} ${NT2_MKL_32}) 
    endif()

    nt2_module_tool(is_multicore RESULT_VARIABLE RESULT_VAR OUTPUT_QUIET)
    if(RESULT_VAR EQUAL -1 OR RESULT_VAR EQUAL 1)
      set(NT2_ARCH_MULTICORE FALSE)
      find_library(NT2_MKL_SEQ NAMES mkl_sequential mkl_sequential_dll
                   PATHS /opt/intel/mkl/lib/ ${NT2_BLAS_ROOT}
                   PATH_SUFFIXES ia32 intel64
                   )
      set(NT2_BLAS_LIBRARIES ${NT2_BLAS_LIBRARIES} ${NT2_MKL_SEQ})
    else()
      if(NT2_COMPILER_MSVC)
        find_library(NT2_MKL_INTEL_THREAD NAMES mkl_intel_thread_dll
                     PATHS /opt/intel/mkl/lib/ ${NT2_BLAS_ROOT}
                     PATH_SUFFIXES ia32 intel64
                     )
        set(NT2_BLAS_LIBRARIES ${NT2_BLAS_LIBRARIES} ${NT2_ICC_LIB_ROOT}/libiomp5md.lib)
      elseif(NT2_COMPILER_GCC)
        find_library(NT2_MKL_GNU_THREAD NAMES mkl_gnu_thread
                     PATHS /opt/intel/mkl/lib/ ${NT2_BLAS_ROOT}
                     PATH_SUFFIXES ia32 intel64
                    ) 
        set(NT2_BLAS_LIBRARIES ${NT2_BLAS_LIBRARIES} ${NT2_MKL_GNU_THREAD})
        set(NT2_BLAS_LINK_FLAGS ${NT2_BLAS_LINK_FLAGS} ${NT2_OPENMP_LINK_FLAGS})
      elseif(NT2_COMPILER_ICC)
        find_library(NT2_MKL_INTEL_THREAD NAMES mkl_intel_thread mkl_intel_thread_dll
                     PATHS /opt/intel/mkl/lib/ ${NT2_BLAS_ROOT}
                     PATH_SUFFIXES ia32 intel64
                    )
        set(NT2_BLAS_LIBRARIES ${NT2_BLAS_LIBRARIES} ${NT2_MKL_INTEL_THREAD})
        if(UNIX)
          set(NT2_BLAS_LINK_FLAGS ${NT2_BLAS_LINK_FLAGS} "-openmp")
        elseif(WIN32)
          set(NT2_BLAS_LINK_FLAGS ${NT2_BLAS_LINK_FLAGS} "/Qopenmp")
        endif()
        set(NT2_ARCH_MULTICORE TRUE)      
      endif()
    endif()

    
    
    if(UNIX)
      find_package(Threads)
      set(NT2_BLAS_LIBRARIES ${NT2_BLAS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
      set(NT2_BLAS_LIBRARIES ${NT2_BLAS_LIBRARIES} "-lm")
    endif()

    find_library(NT2_MKL_CORE NAMES mkl_core mkl_core_dll
                 PATHS /opt/intel/mkl/lib/ ${NT2_BLAS_ROOT}
                 PATH_SUFFIXES ia32 intel64
                )
    set(NT2_BLAS_LIBRARIES ${NT2_BLAS_LIBRARIES} ${NT2_MKL_CORE}) 
    
    if(NT2_MKL_LP64 OR NT2_MKL_32)
      if(NT2_MKL_SEQ AND NT2_MKL_CORE)
        set(NT2_BLAS_FOUND TRUE)
      elseif(NT2_MKL_GNU_THREAD OR NT2_MKL_INTEL_THREAD)
        if(NT2_MKL_CORE)
          set(NT2_BLAS_FOUND TRUE)
        endif()
      endif()
    endif()

  elseif(NT2_BLAS_VENDOR STREQUAL "Goto2")
    find_library(NT2_GOTO2 NAMES goto2
                 PATHS ${NT2_BLAS_ROOT})
    if(NT2_GOTO2)
      set(NT2_BLAS_FOUND TRUE)
      set(NT2_BLAS_LIBRARIES ${NT2_BLAS_LIBRARIES} ${NT2_GOTO2})
    endif()

  elseif(NT2_BLAS_VENDOR STREQUAL "Netlib")  
    find_library(NT2_NETLIB NAMES blas_LINUX 
                 PATHS ${NT2_BLAS_ROOT})
    if(NT2_NETLIB)
      set(NT2_BLAS_FOUND TRUE)
      set(NT2_BLAS_LIBRARIES ${NT2_BLAS_LIBRARIES} ${NT2_NETLIB})
    endif()

  else() # TODO : Other support if findBlas failed

    if(BLAS_FOUND)
      find_package(BLAS QUIET)
        set(NT2_BLAS_FOUND TRUE)
        set(NT2_BLAS_LIBRARIES ${BLAS_LIBRARIES})
        set(NT2_BLAS_LINK_FLAGS ${BLAS_LINKER_FLAGS})
      
    endif()
  endif()

  if(NOT NT2_BLAS_FOUND)
    message("-- Warning: If Blas is available on the target,")
    message("-- please specify the location and the vendor")
    message("-- by defining NT2_BLAS_ROOT and NT2_BLAS_VENDOR.")
    message("-- (SUPPORTED : \"Intel\",\"Goto2\",\"Netlib\")")
  endif()
