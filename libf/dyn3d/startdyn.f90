MODULE startdyn

  ! From startvar.F, version 1.4
  ! 2006/01/27 15:14:22 Fairhead

  IMPLICIT NONE

  private
  public start_init_dyn, start_inter_3d

  INTEGER fid_dyn, iml_dyn, jml_dyn, llm_dyn, ttm_dyn

  REAL, ALLOCATABLE:: lon_ini(:), lat_ini(:)
  ! (longitude and latitude from the input file, converted to rad)

  real, ALLOCATABLE:: levdyn_ini(:)

CONTAINS

  SUBROUTINE start_init_dyn(tsol_2d, psol)

    ! Host associated variables appearing and modified in this procedure :
    ! iml_dyn, jml_dyn, llm_dyn, ttm_dyn, fid_dyn, lon_ini, lat_ini, levdyn_ini

    USE flincom, only: flininfo, flinopen_nozoom
    use comgeom, only: aire_2d, apoln, apols
    use conf_dat2d_m, only: conf_dat2d
    use inter_barxy_m, only: inter_barxy
    use comgeom, only: rlonu, rlatv
    use dimens_m, only: iim, jjm
    use gr_int_dyn_m, only: gr_int_dyn
    use start_init_orog_m, only: phis
    use start_init_phys_m, only: start_init_phys
    use nr_util, only: assert, pi
    use netcdf, only: nf90_nowrite
    use netcdf95, only: nf95_open, nf95_close, nf95_get_var, nf95_inq_varid

    REAL, intent(in):: tsol_2d(:, :) ! (iim + 1, jjm + 1)
    REAL, intent(out):: psol(:, :) ! (iim + 1, jjm + 1) surface pressure, in Pa

    ! Local:

    REAL date, dt
    INTEGER itau(1), ncid, varid
    REAL, ALLOCATABLE:: lon_rad(:), lat_rad(:)

    REAL, ALLOCATABLE:: lon_dyn(:, :), lat_dyn(:, :)
    ! (longitude and latitude from the input file, in rad or degrees)

    REAL, ALLOCATABLE:: var_ana(:, :)
    real z(iim + 1, jjm + 1)
    real tmp_var(iim, jjm + 1)

    !--------------------------

    print *, "Call sequence information: start_init_dyn"
    call assert((/size(tsol_2d, 1), size(psol, 1)/) == iim + 1, &
         "start_init_dyn size 1")
    call assert((/size(tsol_2d, 2), size(psol, 2)/) == jjm + 1, &
         "start_init_dyn size 2")
    CALL flininfo('ECDYN.nc', iml_dyn, jml_dyn, llm_dyn, ttm_dyn, fid_dyn)
    print *, "iml_dyn = ", iml_dyn, ", jml_dyn = ", jml_dyn, &
         ", llm_dyn = ", llm_dyn, ", ttm_dyn = ", ttm_dyn

    ALLOCATE(lat_dyn(iml_dyn, jml_dyn))
    ALLOCATE(lon_dyn(iml_dyn, jml_dyn))
    ALLOCATE(levdyn_ini(llm_dyn))

    CALL flinopen_nozoom(iml_dyn, jml_dyn, llm_dyn, &
         lon_dyn, lat_dyn, levdyn_ini, ttm_dyn, itau, date, dt, fid_dyn)

    ALLOCATE(var_ana(iml_dyn, jml_dyn))
    ALLOCATE(lon_rad(iml_dyn))
    ALLOCATE(lon_ini(iml_dyn))

    IF (MAXVAL(lon_dyn) > pi) THEN
       ! Assume "lon_dyn" is in degrees 
       lon_ini = lon_dyn(:, 1) * pi / 180.
    ELSE
       lon_ini = lon_dyn(:, 1) 
    ENDIF

    ALLOCATE(lat_rad(jml_dyn))
    ALLOCATE(lat_ini(jml_dyn))

    IF (MAXVAL(lat_dyn) > pi) THEN
       lat_ini = lat_dyn(1, :) * pi / 180.
    ELSE
       lat_ini = lat_dyn(1, :) 
    ENDIF

    call nf95_open('ECDYN.nc', nf90_nowrite, ncid)

    ! 'Z': Surface geopotential
    call nf95_inq_varid(ncid, 'Z', varid)
    call nf95_get_var(ncid, varid, var_ana)
    CALL conf_dat2d(lon_ini, lat_ini, lon_rad, lat_rad, var_ana)
    CALL inter_barxy(lon_rad, lat_rad(:jml_dyn -1), var_ana, rlonu(:iim), &
         rlatv, tmp_var)
    z = gr_int_dyn(tmp_var)

    ! 'SP': Surface pressure
    call nf95_inq_varid(ncid, 'SP', varid)
    call nf95_get_var(ncid, varid, var_ana)
    CALL conf_dat2d(lon_ini, lat_ini, lon_rad, lat_rad, var_ana)
    CALL inter_barxy(lon_rad, lat_rad(:jml_dyn -1), var_ana, rlonu(:iim), &
         rlatv, tmp_var) 
    psol = gr_int_dyn(tmp_var)

    call nf95_close(ncid)

    psol(:iim, :) = psol(:iim, :) &
         * (1. + (z(:iim, :) - phis(:iim, :)) / 287. / tsol_2d(:iim, :))
    psol(iim + 1, :) = psol(1, :)

    psol(:, 1) = SUM(aire_2d(:iim, 1) * psol(:iim, 1)) / apoln
    psol(:, jjm + 1) = SUM(aire_2d(:iim, jjm + 1) * psol(:iim, jjm + 1)) &
         / apols

  END SUBROUTINE start_init_dyn

  !********************************

  subroutine start_inter_3d(varname, lon_in2, lat_in2, pls_in, var3d)

    ! This procedure gets a 3D variable from a file and interpolates it.

    use nr_util, only: assert_eq
    use numer_rec, only: spline, splint
    use inter_barxy_m, only: inter_barxy
    use gr_int_dyn_m, only: gr_int_dyn
    use conf_dat3d_m, only: conf_dat3d
    use netcdf, only: nf90_nowrite
    use netcdf95, only: nf95_open, nf95_close, nf95_get_var, nf95_inq_varid

    CHARACTER(len=*), intent(in):: varname
    REAL, intent(in):: lon_in2(:) ! (iml)
    REAL, intent(in):: lat_in2(:)
    REAL, intent(in):: pls_in(:, :, :) ! (iml, jml, lml)
    REAL, intent(out):: var3d(:, :, :) ! (iml, jml, lml)

    ! LOCAL:
    INTEGER iml, jml, lml, ncid, varid
    INTEGER ii, ij, il
    REAL lon_rad(iml_dyn), lat_rad(jml_dyn)
    REAL lev_dyn(llm_dyn)
    REAL var_tmp2d(size(lon_in2)-1, size(pls_in, 2))
    real var_tmp3d(size(lon_in2), size(pls_in, 2), llm_dyn)
    REAL ax(llm_dyn), ay(llm_dyn), yder(llm_dyn)
    real var_ana3d(iml_dyn, jml_dyn, llm_dyn)

    !--------------------------------

    print *, "Call sequence information: start_inter_3d"

    iml = assert_eq(size(pls_in, 1), size(lon_in2), size(var3d, 1), &
         "start_inter_3d iml")
    jml = assert_eq(size(pls_in, 2), size(var3d, 2), "start_inter_3d jml")
    lml = assert_eq(size(pls_in, 3), size(var3d, 3), "start_inter_3d lml")

    print *, "iml = ", iml, ", jml = ", jml
    print *, "varname = ", varname
    call nf95_open('ECDYN.nc', nf90_nowrite, ncid)
    call nf95_inq_varid(ncid, varname, varid)
    call nf95_get_var(ncid, varid, var_ana3d)
    call nf95_close(ncid)
    CALL conf_dat3d(lon_ini, lat_ini, levdyn_ini, lon_rad, lat_rad, lev_dyn, &
         var_ana3d)

    DO il = 1, llm_dyn
       CALL inter_barxy(lon_rad, lat_rad(:jml_dyn-1), var_ana3d(:, :, il), &
            lon_in2(:iml-1), lat_in2, var_tmp2d) 
       var_tmp3d(:, :, il) = gr_int_dyn(var_tmp2d)
    ENDDO

    ! Pour l'interpolation verticale, on interpole du haut de l'atmosphère
    ! vers le sol :
    ax = lev_dyn(llm_dyn:1:-1)
    DO ij=1, jml
       DO ii=1, iml-1
          ay = var_tmp3d(ii, ij, llm_dyn:1:-1)
          yder = SPLINE(ax, ay)
          do il=1, lml
             var3d(ii, ij, il) = SPLINT(ax, ay, yder, pls_in(ii, ij, il))
          END do
       ENDDO
    ENDDO
    var3d(iml, :, :) = var3d(1, :, :) 

  END subroutine start_inter_3d

END MODULE startdyn
