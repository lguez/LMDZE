program test_inifilr

  use dimensions, only: iim, jjm
  use dynetat0_m, only: rlatu, rlatv, read_serre, fyhyp, fxhyp
  use inifilr_m, only: inifilr, jfiltnu, jfiltnv, jfiltsu, jfiltsv, &
       matriceun, matrinvn, matricevn, matriceus, matrinvs, matricevs
  use netcdf, only: NF90_CLOBBER, NF90_FLOAT
  use netcdf95, only: nf95_create, NF95_DEF_DIM, NF95_DEF_VAR, NF95_ENDDEF, &
       NF95_PUT_VAR, NF95_CLOSE, nf95_put_att
  use nr_util, only: pi
  use unit_nml_m, only: unit_nml, set_unit_nml

  IMPLICIT NONE

  integer ncid, dimid_longitude
  integer dimid_rlatu, dimid_rlatv
  integer varid_matriceu, varid_matrinv, varid_matricev
  integer varid_rlatu, varid_rlatv

  !-----------------------------------------------------------------------

  call set_unit_nml
  open(unit_nml, file="used_namelists.txt", status="replace", action="write")
  call read_serre
  CALL fyhyp
  CALL fxhyp
  CALL inifilr
  close(unit_nml)

  call nf95_create("test_inifilr_out.nc", NF90_CLOBBER, ncid)

  call NF95_DEF_DIM(ncid, "longitude", iim, dimid_longitude)
  call NF95_DEF_DIM(ncid, "rlatu", jfiltnu + jjm - jfiltsu, dimid_rlatu)
  call NF95_DEF_DIM(ncid, "rlatv", jfiltnv + jjm - jfiltsv + 1, dimid_rlatv)

  call NF95_DEF_VAR(ncid, "rlatu", NF90_FLOAT, dimid_rlatu, varid_rlatu)
  call nf95_put_att(ncid, varid_rlatu, "units", "degrees_north")

  call NF95_DEF_VAR(ncid, "rlatv", NF90_FLOAT, dimid_rlatv, varid_rlatv)
  call nf95_put_att(ncid, varid_rlatv, "units", "degrees_north")

  call NF95_DEF_VAR(ncid, "matriceu", NF90_FLOAT, &
       (/dimid_longitude, dimid_longitude, dimid_rlatu/), varid_matriceu)
  call NF95_DEF_VAR(ncid, "matrinv", NF90_FLOAT, &
       (/dimid_longitude, dimid_longitude, dimid_rlatu/), varid_matrinv)
  call NF95_DEF_VAR(ncid, "matricev", NF90_FLOAT, &
       (/dimid_longitude, dimid_longitude, dimid_rlatv/), varid_matricev)

  call NF95_ENDDEF(ncid)

  call NF95_PUT_VAR(ncid, varid_rlatu, &
       (/rlatu(2:jfiltnu), rlatu(jfiltsu:jjm)/) / pi * 180.)
  call NF95_PUT_VAR(ncid, varid_rlatv, &
       (/rlatv(:jfiltnv), rlatv(jfiltsv:)/) / pi * 180.)

  call NF95_PUT_VAR(ncid, varid_matriceu, matriceun)
  call NF95_PUT_VAR(ncid, varid_matrinv, matrinvn)
  call NF95_PUT_VAR(ncid, varid_matricev, matricevn)
  call NF95_PUT_VAR(ncid, varid_matriceu, matriceus, start = (/1, 1, jfiltnu/))
  call NF95_PUT_VAR(ncid, varid_matrinv, matrinvs, start = (/1, 1, jfiltnu/))
  call NF95_PUT_VAR(ncid, varid_matricev, matricevs, &
       start = (/1, 1, jfiltnv + 1/))

  call NF95_CLOSE(ncid)
  print *, 'Created file "test_inifilr_out.nc".'

end program test_inifilr
