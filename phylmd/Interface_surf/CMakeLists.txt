target_sources(LMDZE_objects PRIVATE alboc_cd.f90 albsno.f90
  calcul_fluxs.f90 cdrag.f90 climb_hq_down.f90 climb_hq_up.f90
  clqh.f90 clvent.f90 coef_diff_turb.f90 coefkz2.f90 coefkz.f90
  coefkzmin.f90 conf_interface.f90 fonte_neige.f90 hbtm.f90
  interfoce_lim.f90 interfsur_lim.f90 limit_read_sst.f90
  pbl_surface.f90 screenc.f90 screenp.f90 soil.f90 stdlevvar.f90
  ustarhb.f90 yamada4.f90 screencp.f90 calc_coef.f90 dimsoil.f90)

add_executable(test_soil test_soil.f90)

target_link_libraries(test_soil PRIVATE LMDZE_objects Jumble::jumble)
