cmake_minimum_required(VERSION 3.16)
list(APPEND CMAKE_MESSAGE_CONTEXT LMDZE)
project(LMDZE VERSION 0.16 LANGUAGES Fortran CXX)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
option(FETCH "Download Jumble, NetCDF95 and Numer_Rec_95")

if(FETCH)
  include(FetchContent)
  FetchContent_Declare(Jumble GIT_REPOSITORY
    https://github.com/lguez/Jumble.git)
  FetchContent_Declare(NetCDF95 GIT_REPOSITORY
    https://github.com/lguez/NetCDF95.git)
  FetchContent_Declare(Numer_Rec_95 GIT_REPOSITORY
    https://gitlab.in2p3.fr/guez/Numer_Rec_95.git TLS_VERIFY false
    GIT_CONFIG http.sslVerify=false)
  FetchContent_MakeAvailable(Jumble NetCDF95 Numer_Rec_95)
else()
  find_package(Jumble REQUIRED CONFIG)
  find_package(Numer_Rec_95 REQUIRED CONFIG)
  find_package(NetCDF95 REQUIRED CONFIG)
endif()

find_package(NetCDF_Fortran REQUIRED)
find_package(MPI REQUIRED COMPONENTS Fortran CXX)
find_package(XIOS REQUIRED)
add_library(LMDZE_objects OBJECT dimensions.f90)
target_include_directories(LMDZE_objects INTERFACE ${PROJECT_BINARY_DIR})
target_link_libraries(LMDZE_objects PRIVATE NetCDF95::netcdf95
  XIOS::xios NetCDF_Fortran::netcdff Jumble::jumble
  Numer_Rec_95::numer_rec_95)

foreach(my_dir IN ITEMS dyn3d filtrez IOIPSL misc phylmd dynphy_lonlat)
  add_subdirectory(${my_dir})
endforeach(my_dir)

include(TAGS)
