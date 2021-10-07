find_path(XIOS_INCLUDE_DIR xios.mod)
find_library(XIOS_LIBRARY xios)
include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(XIOS REQUIRED_VARS XIOS_INCLUDE_DIR
  XIOS_LIBRARY)

if (XIOS_FOUND)
  add_library(XIOS::xios STATIC IMPORTED)

  set_target_properties(XIOS::xios PROPERTIES IMPORTED_LOCATION
    ${XIOS_LIBRARY} IMPORTED_LINK_INTERFACE_LANGUAGES "CXX")
  
  target_include_directories(XIOS::xios INTERFACE ${XIOS_INCLUDE_DIR})

  target_link_libraries(XIOS::xios INTERFACE MPI::MPI_CXX
    NetCDF_Fortran::netcdff)

endif()
