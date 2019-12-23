set(file_list_var advn.f90 advnqy.f90 advnx.f90 advnz.f90 advnqx.f90
  advnqz.f90 advny.f90)

list(TRANSFORM file_list_var PREPEND ${CMAKE_CURRENT_LIST_DIR}/)
target_sources(LMDZE_objects PRIVATE ${file_list_var})
