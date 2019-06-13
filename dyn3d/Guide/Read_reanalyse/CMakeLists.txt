set(file_list_var nat2gcm.f90 pres2lev.f90 read_reanalyse.f90 reanalyse2nat.f90)
list(TRANSFORM file_list_var PREPEND ${CMAKE_CURRENT_LIST_DIR}/)
target_sources(LMDZE_objects PRIVATE ${file_list_var})
