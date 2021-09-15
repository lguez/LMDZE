list(APPEND CMAKE_MESSAGE_CONTEXT FindNetCDF_Fortran)

if (NOT TARGET NetCDF_Fortran::NetCDF_Fortran)
  find_package(PkgConfig REQUIRED)

  pkg_check_modules(NetCDF_Fortran REQUIRED IMPORTED_TARGET GLOBAL
    netcdf-fortran)

  pkg_get_variable(netcdf_fortran_pcfiledir netcdf-fortran pcfiledir)
  message(VERBOSE "Location of .pc file: ${netcdf_fortran_pcfiledir}")

  message(VERBOSE
    "NetCDF-Fortran include directories: ${NetCDF_Fortran_INCLUDE_DIRS}")

  message(VERBOSE "NetCDF-Fortran libraries: ${NetCDF_Fortran_LINK_LIBRARIES}")

  if(PKG_CONFIG_VERSION_STRING VERSION_LESS "0.29.2"
      OR PKG_CONFIG_VERSION_STRING VERSION_GREATER "1.1")
    # pkg-config strips system flags out of cflags. They do not appear
    # with prefix -isystem when pkg-config is run with
    # --cflags-only-I. So cmake does not get them. So we have to
    # duplicate the call to pkg_get_variable that is already in
    # pkg_check_modules.
    # Version > 1.1 means we are using pkgconf from pkgconf.org
    pkg_get_variable(pkg_netcdf_fortran_includedir netcdf-fortran includedir)

    target_include_directories(PkgConfig::NetCDF_Fortran INTERFACE
      ${pkg_netcdf_fortran_includedir})

    message(VERBOSE
      "pkg_netcdf_fortran_includedir: ${pkg_netcdf_fortran_includedir}")
  endif()

  target_link_libraries(PkgConfig::NetCDF_Fortran INTERFACE NetCDF::NetCDF)
  add_library(NetCDF_Fortran::NetCDF_Fortran ALIAS PkgConfig::NetCDF_Fortran)
endif()

list(POP_BACK CMAKE_MESSAGE_CONTEXT)
