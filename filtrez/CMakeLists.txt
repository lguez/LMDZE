file(STRINGS ${CMAKE_CURRENT_LIST_DIR}/file_list file_list_var)
list(TRANSFORM file_list_var PREPEND ${CMAKE_CURRENT_LIST_DIR}/)
target_sources(LMDZE_objects PRIVATE ${file_list_var})

add_executable(test_inifilr $<TARGET_OBJECTS:LMDZE_objects>
  ${CMAKE_CURRENT_LIST_DIR}/test_inifilr.f90)

target_include_directories(test_inifilr PRIVATE
  ${netcdf95_INCLUDE_DIR} ${netcdff_INCLUDE_DIR})

target_link_libraries(test_inifilr ${library_list})
