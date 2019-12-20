cmake_minimum_required(VERSION 3.14)
project(LMDZE Fortran)

find_path(netcdf95_INCLUDE_DIR netcdf95.mod)
find_library(netcdf95_LIBRARY netcdf95)

find_path(netcdff_INCLUDE_DIR netcdf.mod)
find_library(netcdff_LIBRARY netcdff)

find_path(numer_rec_95_INCLUDE_DIR numer_rec_95.mod)
find_library(numer_rec_95_LIBRARY numer_rec_95)

find_path(jumble_INCLUDE_DIR jumble.mod)
find_library(jumble_LIBRARY jumble)

find_path(nr_util_INCLUDE_DIR nr_util.mod)
find_library(nr_util_LIBRARY nr_util)

set(library_list ${netcdf95_LIBRARY} ${netcdff_LIBRARY}
  ${NetCDF_LIBRARY} ${numer_rec_95_LIBRARY} ${jumble_LIBRARY}
  ${nr_util_LIBRARY})

add_library(LMDZE_objects OBJECT dimensions.F90)

set(CPP_IIM 16 CACHE STRING "iim, number of longitudes")
set(CPP_JJM 12 CACHE STRING "jjm, number of latitudes")
set(CPP_LLM 11 CACHE STRING "llm, number of vertical levels")

message(STATUS "Set iim to ${CPP_IIM}, jjm to ${CPP_JJM}, llm to ${CPP_LLM}")

target_compile_definitions(LMDZE_objects PRIVATE CPP_IIM=${CPP_IIM}
  CPP_JJM=${CPP_JJM} CPP_LLM=${CPP_LLM})

target_include_directories(LMDZE_objects PRIVATE
  ${netcdff_INCLUDE_DIR} ${netcdf95_INCLUDE_DIR} ${jumble_INCLUDE_DIR}
  ${nr_util_INCLUDE_DIR} ${numer_rec_95_INCLUDE_DIR})

foreach(my_dir dyn3d filtrez IOIPSL misc phylmd)
  include(${my_dir}/CMakeLists.txt)
endforeach(my_dir)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
include(ConfigureCompilerFlags)
include(TAGS)