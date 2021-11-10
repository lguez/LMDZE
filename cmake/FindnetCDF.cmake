# Adapted from
# https://gitlab.kitware.com/vtk/vtk/-/blob/master/CMake/FindnetCDF.cmake

#[==[
Provides the following variables:

* `netCDF_FOUND`: Whether netCDF was found or not.
* `netCDF_INCLUDE_DIRS`: Include directories necessary to use netCDF.
* `netCDF_LIBRARIES`: Libraries necessary to use netCDF.
* `netCDF_VERSION`: The version of netCDF found.
* `netCDF::netcdf`: A target to use with `target_link_libraries`.
* `netCDF_HAS_PARALLEL`: Whether or not netCDF was found with parallel IO support.
#]==]

list(APPEND CMAKE_MESSAGE_CONTEXT FindnetCDF)

function(FindnetCDF_get_is_parallel_aware include_dir)
  file(STRINGS "${include_dir}/netcdf_meta.h" _netcdf_lines
    REGEX "#define[ \t]+NC_HAS_PARALLEL[ \t]")
  string(REGEX REPLACE ".*NC_HAS_PARALLEL[ \t]*([0-1]+).*" "\\1"
    _netcdf_has_parallel "${_netcdf_lines}")
  if (_netcdf_has_parallel)
    set(netCDF_HAS_PARALLEL TRUE PARENT_SCOPE)
  else()
    set(netCDF_HAS_PARALLEL FALSE PARENT_SCOPE)
  endif()
endfunction()

find_package(PkgConfig QUIET)

if (PkgConfig_FOUND)
  pkg_check_modules(_netCDF QUIET netcdf IMPORTED_TARGET)
  if (_netCDF_FOUND)
    unset(netCDF_DIR CACHE)
    pkg_get_variable(netcdf_pcfiledir netcdf pcfiledir)
    message(VERBOSE "Location of .pc file: ${netcdf_pcfiledir}")
    message(VERBOSE "netCDF include directories: ${_netCDF_INCLUDE_DIRS}")
    message(VERBOSE "netCDF libraries: ${_netCDF_LINK_LIBRARIES}")
    # Forward the variables in a consistent way.
    set(netCDF_FOUND "${_netCDF_FOUND}")
    set(netCDF_INCLUDE_DIRS "${_netCDF_INCLUDE_DIRS}")
    set(netCDF_LIBRARIES "${_netCDF_LIBRARIES}")
    set(netCDF_VERSION "${_netCDF_VERSION}")

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(netCDF
      REQUIRED_VARS netCDF_LIBRARIES
      # This is not required because system-default include paths
      # are not reported by `FindPkgConfig`, so this might be
      # empty. Assume that if we have a library, the include
      # directories are fine (if any) since PkgConfig reported that
      # the package was found.  netCDF_INCLUDE_DIRS
      VERSION_VAR netCDF_VERSION)

    if (NOT TARGET netCDF::netcdf)
      add_library(netCDF::netcdf INTERFACE IMPORTED)
      set_target_properties(netCDF::netcdf PROPERTIES
        INTERFACE_LINK_LIBRARIES "PkgConfig::_netCDF")
    endif ()

    FindnetCDF_get_is_parallel_aware("${_netCDF_INCLUDEDIR}")
  endif ()
endif ()

if (NOT PkgConfig_FOUND OR NOT _netCDF_FOUND)
  find_path(netCDF_INCLUDE_DIR
    NAMES netcdf.h
    DOC "netcdf include directories")
  ##mark_as_advanced(netCDF_INCLUDE_DIR)

  find_library(netCDF_LIBRARY
    NAMES netcdf
    DOC "netcdf library")
  ##mark_as_advanced(netCDF_LIBRARY)

  if (netCDF_INCLUDE_DIR)
    unset(netCDF_DIR CACHE)
    file(STRINGS "${netCDF_INCLUDE_DIR}/netcdf_meta.h" _netcdf_version_lines
      REGEX "#define[ \t]+NC_VERSION_(MAJOR|MINOR|PATCH|NOTE)")
    string(REGEX REPLACE ".*NC_VERSION_MAJOR *\([0-9]*\).*" "\\1"
      _netcdf_version_major "${_netcdf_version_lines}")
    string(REGEX REPLACE ".*NC_VERSION_MINOR *\([0-9]*\).*" "\\1"
      _netcdf_version_minor "${_netcdf_version_lines}")
    string(REGEX REPLACE ".*NC_VERSION_PATCH *\([0-9]*\).*" "\\1"
      _netcdf_version_patch "${_netcdf_version_lines}")
    string(REGEX REPLACE ".*NC_VERSION_NOTE *\"\([^\"]*\)\".*" "\\1"
      _netcdf_version_note "${_netcdf_version_lines}")
    set(netCDF_VERSION
      "${_netcdf_version_major}.${_netcdf_version_minor}.${_netcdf_version_patch}${_netcdf_version_note}")
    unset(_netcdf_version_major)
    unset(_netcdf_version_minor)
    unset(_netcdf_version_patch)
    unset(_netcdf_version_note)
    unset(_netcdf_version_lines)

    FindnetCDF_get_is_parallel_aware("${netCDF_INCLUDE_DIR}")
  endif ()

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(netCDF
    REQUIRED_VARS netCDF_LIBRARY netCDF_INCLUDE_DIR
    VERSION_VAR netCDF_VERSION)

  if (netCDF_FOUND)
    set(netCDF_INCLUDE_DIRS "${netCDF_INCLUDE_DIR}")
    set(netCDF_LIBRARIES "${netCDF_LIBRARY}")

    if (NOT TARGET netCDF::netcdf)
      add_library(netCDF::netcdf UNKNOWN IMPORTED)
      set_target_properties(netCDF::netcdf PROPERTIES
	IMPORTED_LOCATION "${netCDF_LIBRARY}"
	INTERFACE_INCLUDE_DIRECTORIES "${netCDF_INCLUDE_DIR}")
    endif ()
  endif ()
endif()

list(POP_BACK CMAKE_MESSAGE_CONTEXT)
