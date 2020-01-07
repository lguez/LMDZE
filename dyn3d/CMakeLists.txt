file(STRINGS ${CMAKE_CURRENT_LIST_DIR}/file_list file_list_var)
list(TRANSFORM file_list_var PREPEND ${CMAKE_CURRENT_LIST_DIR}/)
target_sources(LMDZE_objects PRIVATE ${file_list_var})

foreach(my_dir Guide Dissipation Inter_barxy Vlsplt ADVN)
  include(${CMAKE_CURRENT_LIST_DIR}/${my_dir}/CMakeLists.txt)
endforeach(my_dir)

add_executable(gcm ${CMAKE_CURRENT_LIST_DIR}/gcm.f90)
target_link_libraries(gcm LMDZE_objects ${netcdff_LIBRARY} ${NetCDF_LIBRARY})

add_executable(test_fxhyp ${CMAKE_CURRENT_LIST_DIR}/test_fxhyp.f90)
target_link_libraries(test_fxhyp LMDZE_objects ${netcdff_LIBRARY}
  ${NetCDF_LIBRARY})

add_executable(ce0l ${CMAKE_CURRENT_LIST_DIR}/ce0l.f90)
target_link_libraries(ce0l LMDZE_objects ${netcdff_LIBRARY} ${NetCDF_LIBRARY})
