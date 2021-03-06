find_path(XIOS_INCLUDE_DIR xios.mod)
find_library(XIOS_LIBRARY xios)
include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(XIOS REQUIRED_VARS XIOS_INCLUDE_DIR
  XIOS_LIBRARY)

if (XIOS_FOUND)
  add_library(XIOS::XIOS STATIC IMPORTED)

  set_target_properties(XIOS::XIOS PROPERTIES IMPORTED_LOCATION
    ${XIOS_LIBRARY} IMPORTED_LINK_INTERFACE_LANGUAGES "CXX")
  
  target_include_directories(XIOS::XIOS INTERFACE ${XIOS_INCLUDE_DIR})

  target_link_libraries(XIOS::XIOS INTERFACE MPI::MPI_CXX
    NetCDF_Fortran::NetCDF_Fortran)
endif()
