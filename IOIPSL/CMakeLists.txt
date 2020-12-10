file(STRINGS ${CMAKE_CURRENT_LIST_DIR}/file_list file_list_var)
list(TRANSFORM file_list_var PREPEND ${CMAKE_CURRENT_LIST_DIR}/)
target_sources(LMDZE_objects PRIVATE ${file_list_var})

  include(${CMAKE_CURRENT_LIST_DIR}/${my_dir}/CMakeLists.txt)
foreach(my_dir IN ITEMS Mathelp Histcom Calendar)
endforeach(my_dir)
