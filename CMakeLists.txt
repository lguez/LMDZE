cmake_minimum_required(VERSION 3.16)
project(LMDZE LANGUAGES Fortran CXX)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

find_package(NR_util REQUIRED CONFIG)
find_package(Jumble REQUIRED CONFIG)
find_package(Numer_Rec_95 REQUIRED CONFIG)
find_package(NetCDF95 REQUIRED CONFIG)
find_package(NetCDF REQUIRED)
find_package(NetCDF_Fortran REQUIRED)
find_package(MPI REQUIRED COMPONENTS Fortran CXX)

find_path(XIOS_INCLUDE_DIR xios.mod REQUIRED)
find_library(XIOS_LIBRARY xios REQUIRED)
add_library(XIOS STATIC IMPORTED)

set_target_properties(XIOS PROPERTIES IMPORTED_LOCATION
  ${XIOS_LIBRARY} IMPORTED_LINK_INTERFACE_LANGUAGES "CXX")

target_include_directories(XIOS INTERFACE ${XIOS_INCLUDE_DIR})

target_link_libraries(XIOS INTERFACE MPI::MPI_CXX
  NetCDF_Fortran::NetCDF_Fortran)

add_library(LMDZE_objects OBJECT dimensions.f90)
target_include_directories(LMDZE_objects INTERFACE ${PROJECT_BINARY_DIR})

target_link_libraries(LMDZE_objects PRIVATE NetCDF95::netcdf95 XIOS
  NetCDF_Fortran::NetCDF_Fortran jumble numer_rec_95 nr_util)

foreach(my_dir IN ITEMS dyn3d filtrez IOIPSL misc phylmd dynphy_lonlat)
  add_subdirectory(${my_dir})
endforeach(my_dir)

include(ConfigureCompilerFlags)
include(TAGS)
