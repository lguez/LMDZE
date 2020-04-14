set(file_list_var calltherm.f90 ctherm.f90 dqthermcell.f90 thermcell.f90)
list(TRANSFORM file_list_var PREPEND ${CMAKE_CURRENT_LIST_DIR}/)
target_sources(LMDZE_objects PRIVATE ${file_list_var})
