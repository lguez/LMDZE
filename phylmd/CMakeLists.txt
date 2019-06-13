file(STRINGS ${CMAKE_CURRENT_LIST_DIR}/file_list file_list_var)
list(TRANSFORM file_list_var PREPEND ${CMAKE_CURRENT_LIST_DIR}/)
target_sources(LMDZE_objects PRIVATE ${file_list_var})

foreach(my_dir Mobidic Thermcell Orography Radlwsw CV30_routines Conflx
    Interface_surf)
  include(${CMAKE_CURRENT_LIST_DIR}/${my_dir}/CMakeLists.txt)
endforeach(my_dir)

add_executable(test_orbite $<TARGET_OBJECTS:LMDZE_objects>
  ${CMAKE_CURRENT_LIST_DIR}/test_orbite.f90)

target_include_directories(test_orbite PRIVATE ${nr_util_INCLUDE_DIR})
target_link_libraries(test_orbite ${library_list})

add_executable(test_ozonecm $<TARGET_OBJECTS:LMDZE_objects>
  ${CMAKE_CURRENT_LIST_DIR}/test_ozonecm.f90)

target_include_directories(test_ozonecm PRIVATE ${nr_util_INCLUDE_DIR}
  ${netcdff_INCLUDE_DIR})

target_link_libraries(test_ozonecm ${library_list})
