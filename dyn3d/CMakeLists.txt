target_sources(LMDZE_objects PRIVATE abort_gcm.f90 adaptdt.f90
  addfi.f90 advect.f90 advtrac.f90 bernoui.f90 bilan_dyn.f90
  caladvtrac.f90 caldyn0.f90 caldyn.f90 comconst.f90 comgeom.f90
  conf_dat2d.f90 conf_dat3d.f90 conf_gcm.f90 convflu.f90 convmas.f90
  covcont.f90 covnat.f90 dist_sphe.f90 disvert.f90 divergf.f90
  diverg_gam.f90 dteta1.f90 dudv1.f90 dudv2.f90 dynetat0_chosen.f90
  dynetat0.f90 dynredem0.f90 dynredem1.f90 enercin.f90 exner_hyb.f90
  extrapol.f90 flumass.f90 geopot.f90 grad.f90
  grilles_gcm_netcdf_sub.f90 gr_int_dyn_m.f90 groupe.f90 groupeun.f90
  gr_u_scal.f90 gr_v_scal.f90 heavyside.f90 infotrac_init.f90
  init_dynzon.f90 inithist.f90 initial0.f90 integrd.f90
  invert_zoom_x.f90 laplacien.f90 laplacien_gam.f90 laplacien_rot.f90
  laplacien_rotgam.f90 leapfrog.f90 massbar.f90 massbarxy.f90
  massdair.f90 minmax.f90 nxgrad.f90 nxgrad_gam.f90 paramet_m.f90
  principal_cshift.f90 qminimum.f90 q_sat.f90 rotatf.f90
  rotat_nfil.f90 scopy.f90 sortvarc.f90 ssum.f90 tourpot.f90
  unit_nml_m.f90 vitvert.f90 vlspltqs.f90 vlxqs.f90 vlyqs.f90
  writehist.f90)

foreach(my_dir IN ITEMS Guide Dissipation Inter_barxy Vlsplt ADVN)
  add_subdirectory(${my_dir})
endforeach(my_dir)

add_executable(gcm gcm.f90)

target_link_libraries(gcm PRIVATE LMDZE_objects MPI::MPI_Fortran
  NetCDF95::netcdf95 XIOS)

add_executable(test_fxhyp test_fxhyp.f90)
target_link_libraries(test_fxhyp PRIVATE LMDZE_objects)
