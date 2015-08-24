program test_inifilr

  use dimens_m, only: iim, jjm
  use dynetat0_m, only: xprimp025, xprimm025, rlatu1, rlatu2, rlatu, rlatv, &
       yprimu1, yprimu2, rlonu, rlonv, xprimu, xprimv
  use fxhyp_m, only: fxhyp
  use fyhyp_m, only: fyhyp
  use inifilr_m, only: inifilr, jfiltnu, jfiltnv, jfiltsu, jfiltsv, &
       matriceun, matrinvn, matricevn, matriceus, matrinvs, matricevs
  use jumble, only: new_unit
  use netcdf, only: NF90_CLOBBER, NF90_FLOAT
  use netcdf95, only: nf95_create, NF95_DEF_DIM, NF95_DEF_VAR, NF95_ENDDEF, &
       NF95_PUT_VAR, NF95_CLOSE, nf95_put_att

  use nr_util, only: pi
  use unit_nml_m, only: unit_nml
  use read_serre_m, only: read_serre

  IMPLICIT NONE

  integer ncid, dimid_longitude
  integer dimid_rlatu_north, dimid_rlatv_north
  integer dimid_rlatu_south, dimid_rlatv_south
  integer varid_matriceun, varid_matrinvn, varid_matricevn
  integer varid_matriceus, varid_matrinvs, varid_matricevs
  integer varid_rlatu_north, varid_rlatv_north
  integer varid_rlatu_south, varid_rlatv_south

  !-----------------------------------------------------------------------

  call new_unit(unit_nml)
  open(unit_nml, file="used_namelists.txt", status="replace", action="write")
  call read_serre

  CALL fyhyp(rlatu, rlatv, rlatu2, yprimu2, rlatu1, yprimu1)
  CALL fxhyp(xprimm025, rlonv, xprimv, rlonu, xprimu, xprimp025)
  rlatu(1) = pi / 2.
  rlatu(jjm + 1) = -rlatu(1)

  CALL inifilr
  close(unit_nml)

  call nf95_create("test_inifilr_out.nc", NF90_CLOBBER, ncid)

  call NF95_DEF_DIM(ncid, "longitude", iim, dimid_longitude)
  call NF95_DEF_DIM(ncid, "rlatu_north", jfiltnu - 1, dimid_rlatu_north)
  call NF95_DEF_DIM(ncid, "rlatv_north", jfiltnv, dimid_rlatv_north)
  call NF95_DEF_DIM(ncid, "rlatu_south", jjm - jfiltsu + 1, dimid_rlatu_south)
  call NF95_DEF_DIM(ncid, "rlatv_south", jjm - jfiltsv + 1, dimid_rlatv_south)

  call NF95_DEF_VAR(ncid, "rlatu_north", NF90_FLOAT, dimid_rlatu_north, &
       varid_rlatu_north)
  call nf95_put_att(ncid, varid_rlatu_north, "units", "degrees_north")

  call NF95_DEF_VAR(ncid, "rlatv_north", NF90_FLOAT, dimid_rlatv_north, &
       varid_rlatv_north)
  call nf95_put_att(ncid, varid_rlatv_north, "units", "degrees_north")

  call NF95_DEF_VAR(ncid, "rlatu_south", NF90_FLOAT, dimid_rlatu_south, &
       varid_rlatu_south)
  call nf95_put_att(ncid, varid_rlatu_south, "units", "degrees_north")

  call NF95_DEF_VAR(ncid, "rlatv_south", NF90_FLOAT, dimid_rlatv_south, &
       varid_rlatv_south)
  call nf95_put_att(ncid, varid_rlatv_south, "units", "degrees_north")


  call NF95_DEF_VAR(ncid, "matriceun", NF90_FLOAT, &
       (/dimid_longitude, dimid_longitude, dimid_rlatu_north/), varid_matriceun)
  call NF95_DEF_VAR(ncid, "matrinvn", NF90_FLOAT, &
       (/dimid_longitude, dimid_longitude, dimid_rlatu_north/), varid_matrinvn)
  call NF95_DEF_VAR(ncid, "matricevn", NF90_FLOAT, &
       (/dimid_longitude, dimid_longitude, dimid_rlatv_north/), varid_matricevn)
  call NF95_DEF_VAR(ncid, "matriceus", NF90_FLOAT, &
       (/dimid_longitude, dimid_longitude, dimid_rlatu_south/), varid_matriceus)
  call NF95_DEF_VAR(ncid, "matrinvs", NF90_FLOAT, &
       (/dimid_longitude, dimid_longitude, dimid_rlatu_south/), varid_matrinvs)
  call NF95_DEF_VAR(ncid, "matricevs", NF90_FLOAT, &
       (/dimid_longitude, dimid_longitude, dimid_rlatv_south/), varid_matricevs)

  call NF95_ENDDEF(ncid)

  call NF95_PUT_VAR(ncid, varid_rlatu_north, rlatu(2:jfiltnu) / pi * 180.)
  call NF95_PUT_VAR(ncid, varid_rlatv_north, rlatv(:jfiltnv) / pi * 180.)
  call NF95_PUT_VAR(ncid, varid_rlatu_south, rlatu(jfiltsu:jjm) / pi * 180.)
  call NF95_PUT_VAR(ncid, varid_rlatv_south, rlatv(jfiltsv:) / pi * 180.)

  call NF95_PUT_VAR(ncid, varid_matriceun, matriceun)
  call NF95_PUT_VAR(ncid, varid_matrinvn, matrinvn)
  call NF95_PUT_VAR(ncid, varid_matricevn, matricevn)
  call NF95_PUT_VAR(ncid, varid_matriceus, matriceus)
  call NF95_PUT_VAR(ncid, varid_matrinvs, matrinvs)
  call NF95_PUT_VAR(ncid, varid_matricevs, matricevs)

  call NF95_CLOSE(ncid)

end program test_inifilr
