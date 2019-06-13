set(file_list_var vlsplt.f90 vlx.f90 vly.f90 vlz.f90)
list(TRANSFORM file_list_var PREPEND ${CMAKE_CURRENT_LIST_DIR}/)
target_sources(LMDZE_objects PRIVATE ${file_list_var})
