target_sources(LMDZE_objects PRIVATE calfis.f90 etat0phys.f90
  grid_change.f90 grille_m.f90 limit.f90 start_init_dyn.f90
  start_init_phys.f90 start_inter_3d.f90 test_disvert.f90
  etat0dyn.f90 start_init_subsurf.f90)

add_executable(ce0l ce0l.f90)
target_link_libraries(ce0l PRIVATE LMDZE_objects)
