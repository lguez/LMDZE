file(STRINGS ${CMAKE_CURRENT_LIST_DIR}/file_list file_list_var)
list(TRANSFORM file_list_var PREPEND ${CMAKE_CURRENT_LIST_DIR}/)
target_sources(LMDZE_objects PRIVATE ${file_list_var})

add_executable(test_inter_barxy ${CMAKE_CURRENT_LIST_DIR}/test_inter_barxy.f90)
target_link_libraries(test_inter_barxy LMDZE_objects ${netcdff_LIBRARY}
  ${NetCDF_LIBRARY})
