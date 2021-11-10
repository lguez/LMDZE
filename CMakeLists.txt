cmake_minimum_required(VERSION 3.16)
list(APPEND CMAKE_MESSAGE_CONTEXT LMDZE)
project(LMDZE LANGUAGES Fortran CXX)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
option(FETCH "Download")

if(FETCH)
  include(FetchContent)
  
  FetchContent_Declare(NR_util GIT_REPOSITORY
    https://github.com/lguez/NR_util.git)

  FetchContent_Declare(Jumble GIT_REPOSITORY
    https://github.com/lguez/Jumble.git)

  FetchContent_Declare(NetCDF95 GIT_REPOSITORY
    https://github.com/lguez/NetCDF95.git)

  FetchContent_Declare(Numer_Rec_95 GIT_REPOSITORY
    https://gitlab.in2p3.fr/guez/Numer_Rec_95.git)

  FetchContent_MakeAvailable(NR_util Jumble NetCDF95 Numer_Rec_95)
else()
  find_package(NR_util REQUIRED CONFIG)
  find_package(Jumble REQUIRED CONFIG)
  find_package(Numer_Rec_95 REQUIRED CONFIG)
  find_package(NetCDF95 REQUIRED CONFIG)
  find_package(NetCDF REQUIRED)
  find_package(NetCDF_Fortran REQUIRED)
endif()

find_package(MPI REQUIRED COMPONENTS Fortran CXX)
find_package(XIOS REQUIRED)

add_library(LMDZE_objects OBJECT dimensions.f90)
target_include_directories(LMDZE_objects INTERFACE ${PROJECT_BINARY_DIR})

target_link_libraries(LMDZE_objects PRIVATE NetCDF95::netcdf95
  XIOS::xios NetCDF_Fortran::netcdff Jumble::jumble
  Numer_Rec_95::numer_rec_95 NR_util::nr_util)

foreach(my_dir IN ITEMS dyn3d filtrez IOIPSL misc phylmd dynphy_lonlat)
  add_subdirectory(${my_dir})
endforeach(my_dir)

include(TAGS)
