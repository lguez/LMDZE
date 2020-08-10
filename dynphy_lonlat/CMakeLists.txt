set(file_list_var calfis.f90 etat0.f90 grid_change.f90 grille_m.f90
  limit.f90 start_init_dyn.f90 start_init_phys_m.f90
  start_inter_3d.f90 test_disvert.f90)

list(TRANSFORM file_list_var PREPEND ${CMAKE_CURRENT_LIST_DIR}/)
target_sources(LMDZE_objects PRIVATE ${file_list_var})

add_executable(ce0l ${CMAKE_CURRENT_LIST_DIR}/ce0l.f90)
target_link_libraries(ce0l LMDZE_objects)
