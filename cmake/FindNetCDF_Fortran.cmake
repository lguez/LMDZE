list(APPEND CMAKE_MESSAGE_CONTEXT FindNetCDF_Fortran)

if(TARGET NetCDF_Fortran::netcdff)
  set(NetCDF_Fortran_FOUND True)
else()
  # Find NetCDF dependency:

  option(use_find_netcdf_module "Use the find module for NetCDF")

  unset(extraArgs)

  if(${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
    list(APPEND extraArgs QUIET)
  endif()
  
  if(${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED)
    list(APPEND extraArgs REQUIRED)
  endif()
  
  if(use_find_netcdf_module)
    find_package(netCDF ${extraArgs})
  else()
    find_package(netCDF CONFIG ${extraArgs})
  endif()

  #-

  if(netCDF_FOUND)
    find_package(PkgConfig REQUIRED)

    pkg_check_modules(netcdff REQUIRED IMPORTED_TARGET GLOBAL
      netcdf-fortran)

    pkg_get_variable(netcdf_fortran_pcfiledir netcdf-fortran pcfiledir)
    message(VERBOSE "Location of .pc file: ${netcdf_fortran_pcfiledir}")

    message(VERBOSE
      "NetCDF-Fortran include directories: ${netcdff_INCLUDE_DIRS}")

    message(VERBOSE "NetCDF-Fortran libraries: ${netcdff_LINK_LIBRARIES}")

    if(PKG_CONFIG_VERSION_STRING VERSION_LESS "0.29.2"
	OR PKG_CONFIG_VERSION_STRING VERSION_GREATER "1.1")
      # pkg-config strips system flags out of cflags. They do not appear
      # with prefix -isystem when pkg-config is run with
      # --cflags-only-I. So cmake does not get them. So we have to
      # duplicate the call to pkg_get_variable that is already in
      # pkg_check_modules.
      # Version > 1.1 means we are using pkgconf from pkgconf.org
      pkg_get_variable(pkg_netcdf_fortran_includedir netcdf-fortran includedir)

      target_include_directories(PkgConfig::netcdff INTERFACE
	${pkg_netcdf_fortran_includedir})

      message(VERBOSE
	"pkg_netcdf_fortran_includedir: ${pkg_netcdf_fortran_includedir}")
    endif()

    target_link_libraries(PkgConfig::netcdff INTERFACE netCDF::netcdf)
    add_library(NetCDF_Fortran::netcdff ALIAS PkgConfig::netcdff)
    set(NetCDF_Fortran_FOUND True)
  else()
    set(NetCDF_Fortran_FOUND False)
  endif()
endif()

list(POP_BACK CMAKE_MESSAGE_CONTEXT)
