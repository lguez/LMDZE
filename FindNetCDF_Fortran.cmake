find_package(PkgConfig REQUIRED)
pkg_check_modules(NetCDF_Fortran REQUIRED IMPORTED_TARGET netcdf-fortran)
pkg_get_variable(netcdf_fortran_pcfiledir netcdf-fortran pcfiledir)
message(STATUS "${CMAKE_CURRENT_LIST_FILE}")
message(STATUS "-- Location of .pc file: ${netcdf_fortran_pcfiledir}")

message(STATUS
  "-- NetCDF-Fortran include directories: ${NetCDF_Fortran_INCLUDE_DIRS}")

message(STATUS "-- NetCDF-Fortran libraries: ${NetCDF_Fortran_LINK_LIBRARIES}")

if(PKG_CONFIG_VERSION_STRING VERSION_LESS "0.29.2")
  # pkg-config strips system flags out of cflags. They do not appear
  # with prefix -isystem when pkg-config is run with
  # --cflags-only-I. So cmake does not get them. So we have to
  # duplicate the call to pkg_get_variable that is already in
  # pkg_check_modules.
  pkg_get_variable(pkg_netcdf_fortran_includedir netcdf-fortran includedir)

  target_include_directories(PkgConfig::NetCDF_Fortran INTERFACE
    ${pkg_netcdf_fortran_includedir})

  message(STATUS
    "-- pkg_netcdf_fortran_includedir: ${pkg_netcdf_fortran_includedir}")
endif()

add_library(NetCDF_Fortran INTERFACE IMPORTED)

set_target_properties(NetCDF_Fortran PROPERTIES
  INTERFACE_LINK_LIBRARIES PkgConfig::NetCDF_Fortran)

target_link_libraries(NetCDF_Fortran INTERFACE NetCDF::NetCDF)
add_library(NetCDF_Fortran::NetCDF_Fortran ALIAS NetCDF_Fortran)
