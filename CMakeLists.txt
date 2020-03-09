cmake_minimum_required(VERSION 3.16)
project(LMDZE LANGUAGES Fortran CXX)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})

find_path(netcdff_INCLUDE_DIR netcdf.mod)

find_path(XIOS_INCLUDE_DIR xios.mod)
find_library(XIOS_LIBRARY xios)

find_package(NR_util REQUIRED CONFIG)
find_package(Jumble REQUIRED CONFIG)
find_package(Numer_Rec_95 REQUIRED CONFIG)
find_package(NetCDF95 REQUIRED CONFIG)
find_package(NetCDF REQUIRED)
find_package(MPI REQUIRED COMPONENTS CXX)

add_library(LMDZE_objects OBJECT dimensions.F90)

set(CPP_IIM 16 CACHE STRING "iim, number of longitudes")
set(CPP_JJM 12 CACHE STRING "jjm, number of latitudes")
set(CPP_LLM 11 CACHE STRING "llm, number of vertical levels")

message(STATUS "Set iim to ${CPP_IIM}, jjm to ${CPP_JJM}, llm to ${CPP_LLM}")

target_compile_definitions(LMDZE_objects PRIVATE CPP_IIM=${CPP_IIM}
  CPP_JJM=${CPP_JJM} CPP_LLM=${CPP_LLM})

target_include_directories(LMDZE_objects PRIVATE ${netcdff_INCLUDE_DIR})
target_link_libraries(LMDZE_objects netcdf95 jumble numer_rec_95)

foreach(my_dir dyn3d filtrez IOIPSL misc phylmd)
  include(${my_dir}/CMakeLists.txt)
endforeach(my_dir)

include(ConfigureCompilerFlags)
include(TAGS)
