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

    USE ioipsl, only: flininfo, flinopen_nozoom, flinget
    use comgeom, only: aire_2d, apoln, apols
    use conf_dat2d_m, only: conf_dat2d
    use inter_barxy_m, only: inter_barxy
    use comconst, only: pi
    use comgeom, only: rlonu, rlatv
    use dimens_m, only: iim, jjm
    use gr_int_dyn_m, only: gr_int_dyn
    use start_init_orog_m, only: phis
    use start_init_phys_m, only: start_init_phys

    REAL, intent(out):: tsol_2d(:, :)
    REAL, intent(out):: psol(:, :) ! surface pressure, in Pa

    ! Local:

    REAL date, dt
    INTEGER itau(1)
    INTEGER i, j

    CHARACTER(len=120) physfname

    REAL, ALLOCATABLE:: lon_rad(:), lat_rad(:)

    REAL, ALLOCATABLE:: lon_dyn(:, :), lat_dyn(:, :)
    ! (longitude and latitude from the input file, in rad or degrees)

    REAL, ALLOCATABLE:: var_ana(:, :), z(:, :)
    real tmp_var(iim, jjm + 1)
    REAL, ALLOCATABLE:: xppn(:), xpps(:)

    !--------------------------

    print *, "Call sequence information: start_init_dyn"
    if (any((/size(tsol_2d, 1), size(psol, 1)/) /= iim + 1)) stop &
         "start_init_phys size 1"
    if (any((/size(tsol_2d, 2), size(psol, 2)/) /= jjm + 1)) stop &
         "start_init_phys size 2"
    physfname = 'ECDYN.nc'
    print *, 'Opening the surface analysis'
    CALL flininfo(physfname, iml_dyn, jml_dyn, llm_dyn, ttm_dyn, fid_dyn)
    print *, 'Values read from "' // trim(physfname) // '":'
    print *, "iml_dyn = ", iml_dyn, ", jml_dyn = ", jml_dyn, &
         ", llm_dyn = ", llm_dyn, ", ttm_dyn = ", ttm_dyn

    ALLOCATE(lat_dyn(iml_dyn, jml_dyn))
    ALLOCATE(lon_dyn(iml_dyn, jml_dyn))
    ALLOCATE(levdyn_ini(llm_dyn))

    CALL flinopen_nozoom(physfname, iml_dyn, jml_dyn, llm_dyn, &
         lon_dyn, lat_dyn, levdyn_ini, ttm_dyn, itau, date, dt, fid_dyn)

    ALLOCATE(var_ana(iml_dyn, jml_dyn))
    ALLOCATE(lon_rad(iml_dyn))
    ALLOCATE(lon_ini(iml_dyn))

    IF (MAXVAL(lon_dyn(:, :)) > pi) THEN
       ! Assume "lon_dyn" is in degrees 
       lon_ini(:) = lon_dyn(:, 1) * pi / 180.
    ELSE
       lon_ini(:) = lon_dyn(:, 1) 
    ENDIF

    ALLOCATE(lat_rad(jml_dyn))
    ALLOCATE(lat_ini(jml_dyn))

    IF (MAXVAL(lat_dyn(:, :)) > pi) THEN
       lat_ini(:) = lat_dyn(1, :) * pi / 180.
    ELSE
       lat_ini(:) = lat_dyn(1, :) 
    ENDIF

    ALLOCATE(z(iim + 1, jjm + 1))

    ! 'Z': Surface geopotential
    CALL flinget(fid_dyn, 'Z', iml_dyn, jml_dyn, 0, ttm_dyn, 1, 1, var_ana)
    CALL conf_dat2d(lon_ini, lat_ini, lon_rad, lat_rad, var_ana)
    CALL inter_barxy(lon_rad, lat_rad(:jml_dyn -1), var_ana, rlonu(:iim), &
         rlatv, tmp_var)
    z(:, :) = gr_int_dyn(tmp_var)

    ! 'SP': Surface pressure
    CALL flinget(fid_dyn, 'SP', iml_dyn, jml_dyn, 0, ttm_dyn, 1, 1, var_ana)
    CALL conf_dat2d(lon_ini, lat_ini, lon_rad, lat_rad, var_ana)
    CALL inter_barxy(lon_rad, lat_rad(:jml_dyn -1), var_ana, rlonu(:iim), &
         rlatv, tmp_var) 
    psol(:, :) = gr_int_dyn(tmp_var)
    CALL start_init_phys(tsol_2d)

    ! PSOL is computed in Pascals

    DO j = 1, jjm + 1
       DO i = 1, iim
          psol(i, j) = psol(i, j) &
               * (1. + (z(i, j) - phis(i, j)) / 287. / tsol_2d(i, j))
       ENDDO
    ENDDO
    psol(iim + 1, :) = psol(1, :)

    ALLOCATE(xppn(iim))
    ALLOCATE(xpps(iim)) 

    DO  i = 1, iim
       xppn(i) = aire_2d( i, 1) * psol( i, 1)
       xpps(i) = aire_2d( i, jjm + 1) * psol( i, jjm + 1)
    ENDDO

    psol(:, 1) = SUM(xppn)/apoln
    psol(:, jjm + 1) = SUM(xpps)/apols

  END SUBROUTINE start_init_dyn

  !********************************

  function start_inter_3d(varname, lon_in2, lat_in2, pls_in)

    ! This procedure gets a 3D variable from a file and does the
    ! interpolations needed.

    USE ioipsl, only: flinget
    use numer_rec, only: assert_eq, spline, splint
    use inter_barxy_m, only: inter_barxy
    use gr_int_dyn_m, only: gr_int_dyn
    use conf_dat3d_m, only: conf_dat3d

    CHARACTER(len=*), intent(in):: varname
    REAL, intent(in):: lon_in2(:), lat_in2(:)
    REAL, intent(in):: pls_in(:, :, :)

    REAL start_inter_3d(size(lon_in2), size(pls_in, 2), size(pls_in, 3))

    ! LOCAL:
    INTEGER iml, jml, lml
    INTEGER ii, ij, il
    REAL lon_rad(iml_dyn), lat_rad(jml_dyn)
    REAL lev_dyn(llm_dyn)
    REAL var_tmp2d(size(lon_in2)-1, size(pls_in, 2))
    real var_tmp3d(size(lon_in2), size(pls_in, 2), llm_dyn)
    REAL ax(llm_dyn), ay(llm_dyn), yder(llm_dyn)
    real var_ana3d(iml_dyn, jml_dyn, llm_dyn)

    !--------------------------------

    print *, "Call sequence information: start_inter_3d"

    iml = assert_eq(size(pls_in, 1), size(lon_in2), "start_inter_3d iml")
    jml = size(pls_in, 2)
    lml = size(pls_in, 3)

    print *, "iml = ", iml, ", jml = ", jml
    print *, "fid_dyn = ", fid_dyn, ", varname = ", varname
    print *, "iml_dyn = ", iml_dyn, ", jml_dyn = ", jml_dyn, &
         ", llm_dyn = ", llm_dyn, ", ttm_dyn = ", ttm_dyn
    print *, 'Going into flinget to extract the 3D field.'
    CALL flinget(fid_dyn, varname, iml_dyn, jml_dyn, llm_dyn, ttm_dyn, 1, 1, &
         var_ana3d)

    CALL conf_dat3d(lon_ini, lat_ini, levdyn_ini, lon_rad, lat_rad, lev_dyn, &
         var_ana3d)

    DO il = 1, llm_dyn
       CALL inter_barxy(lon_rad, lat_rad(:jml_dyn-1), var_ana3d(:, :, il), &
            lon_in2(:iml-1), lat_in2, var_tmp2d) 
       var_tmp3d(:, :, il) = gr_int_dyn(var_tmp2d)
    ENDDO

    ! Pour l'interpolation verticale, on interpole du haut de l'atmosphère
    ! vers le sol :
    ax(:) = lev_dyn(llm_dyn:1:-1)
    DO ij=1, jml
       DO ii=1, iml-1
          ay(:) = var_tmp3d(ii, ij, llm_dyn:1:-1)
          yder(:) = SPLINE(ax, ay)
          do il=1, lml
             start_inter_3d(ii, ij, il) &
                  = SPLINT(ax, ay, yder, pls_in(ii, ij, il))
          END do
       ENDDO
    ENDDO
    start_inter_3d(iml, :, :) = start_inter_3d(1, :, :) 

  END function start_inter_3d

END MODULE startdyn
