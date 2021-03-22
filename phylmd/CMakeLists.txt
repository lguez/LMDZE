target_sources(LMDZE_objects PRIVATE aaam_bud.f90 ajsec.f90 chem.f90
  clesphys2.f90 clesphys.f90 clouds_gno.f90 cltrac.f90 cltracrn.f90
  comfisrtilp.f90 concvl.f90 conf_phys.f90 cv_driver.f90 cvltr.f90
  diagcld1.f90 dimphy.f90 dimsoil.f90 FCTTRE.f90 fisrtilp.f90
  gr_phy_write.f90 hgardfou.f90 histwrite_phy.f90 indicesol.f90
  ini_histins.f90 initrrnpb.f90 minmaxqfi.f90 newmicro.f90 nflxtr.f90
  orbite.f90 ozonecm.f90 phyetat0.f90 phyredem0.f90 phyredem.f90
  physiq.f90 phytrac.f90 radepsi.f90 radiornpb.f90 suphec.f90
  time_phylmdz.f90 transp.f90 transp_lay.f90 yoethf.f90 YOMCST.f90
  zenang.f90)

foreach(my_dir IN ITEMS Mobidic Thermcell Orography Radlwsw CV30_routines Conflx
    Interface_surf)
  add_subdirectory(${my_dir})
endforeach(my_dir)

add_executable(test_orbite test_orbite.f90)
target_link_libraries(test_orbite PRIVATE LMDZE_objects)

add_executable(test_ozonecm test_ozonecm.f90)

target_link_libraries(test_ozonecm PRIVATE LMDZE_objects nr_util
  NetCDF_Fortran::NetCDF_Fortran)
