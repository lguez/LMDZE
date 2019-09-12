if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
  # Fortran language options:
  string(APPEND CMAKE_Fortran_FLAGS " -std=f2003")

  # Error and warning options:
  string(APPEND CMAKE_Fortran_FLAGS
    " -fmax-errors=1 -pedantic -Wall -Wcharacter-truncation -Wunused-parameter -Wno-conversion -Wimplicit-interface -Wimplicit-procedure -Wno-integer-division")

  # Debugging options:
  set(CMAKE_Fortran_FLAGS_DEBUG
    "-fbacktrace -g -ffpe-trap=invalid,zero,overflow")

  # Code generation options:
  string(APPEND CMAKE_Fortran_FLAGS_DEBUG
    " -fcheck=bounds,do,mem,pointer,recursion -finit-real=nan -O0")

  # Optimization options:
  set(CMAKE_Fortran_FLAGS_RELEASE -O3)

  # Hardware model options:
  string(APPEND CMAKE_Fortran_FLAGS_RELEASE " -mcmodel=medium")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  # Language:
  string(APPEND CMAKE_Fortran_FLAGS
    " -noaltparam -stand f03 -standard-semantics -assume nostd_mod_proc_name")
  string(APPEND CMAKE_Fortran_FLAGS_DEBUG
    " -check bounds,format,output_conversion,pointers,stack,uninit")

  # Data:
  string(APPEND CMAKE_Fortran_FLAGS " -auto -mcmodel=medium")
  string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -init=arrays,minus_huge,snan")
  
  # Compiler diagnostics:
  string(APPEND CMAKE_Fortran_FLAGS
    " -warn declarations,general,stderrors,truncated_source,uncalled,unused,usage -traceback -diag-error-limit 1")
  
  # Optimization:
  string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -O0")
  
  # Floating point:
  string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -fp-stack-check -fpe-all=0")
  
  # Debug:
  string(APPEND CMAKE_Fortran_FLAGS_DEBUG
    " -debug full -debug-parameters all -ftrapuv")
endif()
