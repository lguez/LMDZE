cmake_minimum_required(VERSION 3.16)
project(LMDZE LANGUAGES Fortran CXX)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})

include(FindPkgConfig)
pkg_check_modules(netcdf-fortran IMPORTED_TARGET netcdf-fortran)

find_path(XIOS_INCLUDE_DIR xios.mod REQUIRED)
find_library(XIOS_LIBRARY xios REQUIRED)

find_package(NR_util REQUIRED CONFIG)
find_package(Jumble REQUIRED CONFIG)
find_package(Numer_Rec_95 REQUIRED CONFIG)
find_package(NetCDF95 REQUIRED CONFIG)
find_package(NetCDF REQUIRED)
find_package(MPI REQUIRED COMPONENTS CXX)

add_library(LMDZE_objects OBJECT dimensions.f90)
target_include_directories(LMDZE_objects INTERFACE ${PROJECT_BINARY_DIR})
target_link_libraries(LMDZE_objects PUBLIC netcdf95 jumble numer_rec_95 nr_util)

foreach(my_dir IN ITEMS dyn3d filtrez IOIPSL misc phylmd dynphy_lonlat)
  add_subdirectory(${my_dir})
endforeach(my_dir)

include(ConfigureCompilerFlags)
include(TAGS)
