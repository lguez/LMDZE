module limit_mod

  ! This module is clean: no C preprocessor directive, no include line.

  IMPLICIT none

contains

  SUBROUTINE limit

    ! Authors: L. Fairhead, Z. X. Li, P. Le Van

    ! This subroutine creates files containing boundary conditions.
    ! It uses files with climatological data.
    ! Both grids must be regular.

    use dimens_m, only: iim, jjm
    use comconst, only: daysec, dtvr
    use indicesol, only: epsfra, nbsrf, is_ter, is_oce, is_lic, is_sic
    use dimphy, only: klon, zmasq
    use conf_gcm_m, only: day_step
    use comgeom, only: rlonu, rlatv
    use etat0_mod, only: pctsrf
    use start_init_orog_m, only: masque
    use conf_dat2d_m, only: conf_dat2d
    use inter_barxy_m, only: inter_barxy
    use numer_rec, only: spline, splint
    use grid_change, only: dyn_phy

    use netcdf95, only: handle_err, nf95_get_coord, NF95_CLOSE, NF95_DEF_DIM, &
         nf95_enddef, NF95_CREATE, nf95_inq_dimid, nf95_inquire_dimension, &
         nf95_inq_varid, NF95_OPEN
    use netcdf, only: NF90_CLOBBER, nf90_def_var, NF90_FLOAT, NF90_GET_VAR, &
         NF90_GLOBAL, NF90_NOWRITE, NF90_PUT_ATT, NF90_PUT_VAR, &
         NF90_UNLIMITED

    ! Variables local to the procedure:

    LOGICAL:: extrap = .FALSE.
    ! (extrapolation de données, comme pour les SST lorsque le fichier
    ! ne contient pas uniquement des points océaniques)

    REAL phy_alb(klon, 360)
    REAL phy_sst(klon, 360)
    REAL phy_bil(klon, 360)
    REAL phy_rug(klon, 360)
    REAL phy_ice(klon)

    real pctsrf_t(klon, nbsrf, 360) ! composition of the surface

    ! Pour le champ de départ:
    INTEGER imdep, jmdep, lmdep

    REAL, ALLOCATABLE:: dlon(:), dlat(:)
    REAL, pointer:: dlon_ini(:), dlat_ini(:), timeyear(:)
    REAL, ALLOCATABLE:: champ(:, :)
    REAL, ALLOCATABLE:: work(:, :)

    ! Pour le champ interpolé 3D :
    REAL, allocatable:: champtime(:, :, :)
    REAL champan(iim + 1, jjm + 1, 360)

    ! Pour l'inteprolation verticale :
    REAL, allocatable:: yder(:)

    INTEGER ierr

    INTEGER nid, ndim, ntim
    INTEGER dims(2), debut(2)
    INTEGER id_tim
    INTEGER id_SST, id_BILS, id_RUG, id_ALB
    INTEGER id_FOCE, id_FSIC, id_FTER, id_FLIC

    INTEGER i, j, k, l
    INTEGER ncid, varid, dimid

    REAL, parameter:: tmidmonth(12) = (/(15. + 30. * i, i = 0, 11)/)

    namelist /limit_nml/extrap

    !--------------------

    print *, "Call sequence information: limit"

    print *, "Enter namelist 'limit_nml'."
    read (unit=*, nml=limit_nml)
    write(unit=*, nml=limit_nml)

    ! Initializations:
    dtvr = daysec / real(day_step)
    CALL inigeom

    ! Process rugosity:

    PRINT *, 'Processing rugosity...'
    call NF95_OPEN('Rugos.nc', NF90_NOWRITE, ncid)

    call nf95_get_coord(ncid, "longitude", dlon_ini)
    imdep = size(dlon_ini)

    call nf95_get_coord(ncid, "latitude", dlat_ini)
    jmdep = size(dlat_ini)

    call nf95_get_coord(ncid, "temps", timeyear)
    lmdep = size(timeyear)

    ALLOCATE(champ(imdep, jmdep), champtime(iim, jjm + 1, lmdep))
    allocate(dlon(imdep), dlat(jmdep))
    call NF95_INQ_VARID(ncid, 'RUGOS', varid)

    ! Compute "champtime":
    DO  l = 1, lmdep
       ierr = NF90_GET_VAR(ncid, varid, champ, start=(/1, 1, l/))
       call handle_err("NF90_GET_VAR", ierr)

       CALL conf_dat2d(dlon_ini, dlat_ini, dlon, dlat, champ)
       CALL inter_barxy(dlon, dlat(:jmdep -1), LOG(champ), rlonu(:iim), &
            rlatv, champtime(:, :, l))
       champtime(:, :, l) = EXP(champtime(:, :, l))
       where (nint(masque(:iim, :)) /= 1) champtime(:, :, l) = 0.001
    end do

    call NF95_CLOSE(ncid)

    DEALLOCATE(dlon, dlat, champ, dlon_ini, dlat_ini)
    allocate(yder(lmdep))

    DO j = 1, jjm + 1
       DO i = 1, iim
          yder(:) = SPLINE(timeyear, champtime(i, j, :))
          DO k = 1, 360
             champan(i, j, k) = SPLINT(timeyear, champtime(i, j, :), yder, &
                  real(k-1))
          ENDDO
       ENDDO
    ENDDO

    deallocate(timeyear, champtime, yder)
    champan(iim + 1, :, :) = champan(1, :, :)
    forall (k = 1:360) phy_rug(:, k) = pack(champan(:, :, k), dyn_phy)

    ! Process sea ice:

    PRINT *, 'Processing sea ice...'
    call NF95_OPEN('amipbc_sic_1x1.nc', NF90_NOWRITE, ncid)

    call nf95_get_coord(ncid, "longitude", dlon_ini)
    imdep = size(dlon_ini)

    call nf95_get_coord(ncid, "latitude", dlat_ini)
    jmdep = size(dlat_ini)

    call nf95_inq_dimid(ncid, "time", dimid)
    call NF95_INQuire_DIMension(ncid, dimid, len=lmdep)
    print *, 'lmdep = ', lmdep
    ! PM 28/02/2002 : nouvelle coordonnée temporelle, fichiers AMIP
    ! pas en jours
    ! Ici on suppose qu'on a 12 mois (de 30 jours).
    IF (lmdep /= 12) STOP 'Unknown AMIP file: not 12 months?'

    ALLOCATE(champ(imdep, jmdep), champtime(iim, jjm + 1, lmdep))
    ALLOCATE (dlon(imdep), dlat(jmdep))
    call NF95_INQ_VARID(ncid, 'sicbcs', varid)
    DO l = 1, lmdep
       ierr = NF90_GET_VAR(ncid, varid, champ, start=(/1, 1, l/))
       call handle_err("NF90_GET_VAR", ierr)

       CALL conf_dat2d(dlon_ini, dlat_ini, dlon, dlat, champ)
       CALL inter_barxy (dlon, dlat(:jmdep -1), champ, rlonu(:iim), rlatv, &
            champtime(:, :, l))
    ENDDO

    call NF95_CLOSE(ncid)

    DEALLOCATE(dlon, dlat, champ, dlon_ini, dlat_ini)
    PRINT *, 'Interpolation temporelle'
    allocate(yder(lmdep))

    DO j = 1, jjm + 1
       DO i = 1, iim
          yder(:) = SPLINE(tmidmonth, champtime(i, j, :))
          DO k = 1, 360
             champan(i, j, k) = SPLINT(tmidmonth, champtime(i, j, :), yder, &
                  real(k-1))
          ENDDO
       ENDDO
    ENDDO

    deallocate(champtime, yder)

    ! Convert from percentage to normal fraction and keep sea ice
    ! between 0 and 1:
    champan(:iim, :, :) = max(0., (min(1., (champan(:iim, :, :) / 100.))))
    champan(iim + 1, :, :) = champan(1, :, :)

    DO k = 1, 360
       phy_ice(:) = pack(champan(:, :, k), dyn_phy)

       ! (utilisation de la sous-maille fractionnelle tandis que l'ancien
       ! codage utilisait l'indicateur du sol (0, 1, 2, 3))
       ! PB en attendant de mettre fraction de terre
       WHERE(phy_ice(:) < EPSFRA) phy_ice(:) = 0.

       pctsrf_t(:, is_ter, k) = pctsrf(:, is_ter)
       pctsrf_t(:, is_lic, k) = pctsrf(:, is_lic)
       pctsrf_t(:, is_sic, k) = max(phy_ice(:) - pctsrf_t(:, is_lic, k), 0.)
       ! Il y a des cas où il y a de la glace dans landiceref et
       ! pas dans AMIP
       WHERE( 1. - zmasq(:) < EPSFRA)
          pctsrf_t(:, is_sic, k) = 0.
          pctsrf_t(:, is_oce, k) = 0.
       elsewhere
          where (pctsrf_t(:, is_sic, k) >= 1 - zmasq(:))
             pctsrf_t(:, is_sic, k) = 1. - zmasq(:)
             pctsrf_t(:, is_oce, k) = 0.
          ELSEwhere
             pctsrf_t(:, is_oce, k) = 1. - zmasq(:) - pctsrf_t(:, is_sic, k)
             where (pctsrf_t(:, is_oce, k) < EPSFRA)
                pctsrf_t(:, is_oce, k) = 0.
                pctsrf_t(:, is_sic, k) = 1 - zmasq(:)
             end where
          end where
       end where

       DO i = 1, klon
          if (pctsrf_t(i, is_oce, k) < 0.) then
             print *, 'Problème sous maille : pctsrf_t(', i, &
                  ', is_oce, ', k, ') = ', pctsrf_t(i, is_oce, k)
          ENDIF
          IF (abs(pctsrf_t(i, is_ter, k) + pctsrf_t(i, is_lic, k) &
               + pctsrf_t(i, is_oce, k) + pctsrf_t(i, is_sic, k) - 1.) &
               > EPSFRA) THEN 
             print *, 'Problème sous surface :'
             print *, "pctsrf_t(", i, ", :, ", k, ") = ", &
                  pctsrf_t(i, :, k)
             print *, "phy_ice(", i, ") = ", phy_ice(i)
          ENDIF
       END DO
    ENDDO

    PRINT *, 'Traitement de la sst'
    call NF95_OPEN('amipbc_sst_1x1.nc', NF90_NOWRITE, ncid)

    call nf95_get_coord(ncid, "longitude", dlon_ini)
    imdep = size(dlon_ini)

    call nf95_get_coord(ncid, "latitude", dlat_ini)
    jmdep = size(dlat_ini)

    call nf95_inq_dimid(ncid, "time", dimid)
    call NF95_INQuire_DIMension(ncid, dimid, len=lmdep)
    print *, 'lmdep = ', lmdep
    !PM28/02/2002 : nouvelle coord temporelle fichiers AMIP pas en jours
    !        Ici on suppose qu'on a 12 mois (de 30 jours).
    IF (lmdep /= 12) stop 'Unknown AMIP file: not 12 months?'

    ALLOCATE( champ(imdep, jmdep), champtime(iim, jjm + 1, lmdep))
    IF( extrap )  THEN
       ALLOCATE ( work(imdep, jmdep) )
    ENDIF
    ALLOCATE(   dlon(imdep), dlat(jmdep) )
    call NF95_INQ_VARID(ncid, 'tosbcs', varid)

    DO l = 1, lmdep
       ierr = NF90_GET_VAR(ncid, varid, champ, start=(/1, 1, l/))
       call handle_err("NF90_GET_VAR", ierr)

       CALL conf_dat2d(dlon_ini, dlat_ini, dlon, dlat, champ)
       IF ( extrap ) THEN
          CALL extrapol(champ, imdep, jmdep, 999999., .TRUE., .TRUE., 2, work)
       ENDIF

       CALL inter_barxy (dlon, dlat(:jmdep -1), champ, rlonu(:iim), rlatv, &
            champtime(:, :, l) )
    ENDDO

    call NF95_CLOSE(ncid)

    DEALLOCATE(dlon, dlat, champ, dlon_ini, dlat_ini)
    allocate(yder(lmdep))

    ! interpolation temporelle
    DO j = 1, jjm + 1
       DO i = 1, iim
          yder(:) = SPLINE(tmidmonth, champtime(i, j, :))
          DO k = 1, 360
             champan(i, j, k) = SPLINT(tmidmonth, champtime(i, j, :), yder, &
                  real(k-1))
          ENDDO
       ENDDO
    ENDDO

    deallocate(champtime, yder)
    champan(iim + 1, :, :) = champan(1, :, :)

    !IM14/03/2002 : SST amipbc greater then 271.38
    PRINT *, 'SUB. limit_netcdf.F IM : SST Amipbc >= 271.38 '
    DO k = 1, 360
       DO j = 1, jjm + 1
          DO i = 1, iim
             champan(i, j, k) = amax1(champan(i, j, k), 271.38)
          ENDDO
          champan(iim + 1, j, k) = champan(1, j, k)
       ENDDO
    ENDDO
    forall (k = 1:360) phy_sst(:, k) = pack(champan(:, :, k), dyn_phy)

    ! Traitement de l'albedo

    PRINT *, 'Traitement de l albedo'
    call NF95_OPEN('Albedo.nc', NF90_NOWRITE, ncid)

    call nf95_get_coord(ncid, "longitude", dlon_ini)
    imdep = size(dlon_ini)

    call nf95_get_coord(ncid, "latitude", dlat_ini)
    jmdep = size(dlat_ini)

    call nf95_get_coord(ncid, "temps", timeyear)
    lmdep = size(timeyear)

    ALLOCATE ( champ(imdep, jmdep), champtime(iim, jjm + 1, lmdep))
    ALLOCATE (   dlon(imdep), dlat(jmdep) )
    call NF95_INQ_VARID(ncid, 'ALBEDO', varid)

    DO l = 1, lmdep
       PRINT *, 'Lecture temporelle et int. horizontale ', l, timeyear(l)
       ierr = NF90_GET_VAR(ncid, varid, champ, start=(/1, 1, l/))
       call handle_err("NF90_GET_VAR", ierr)

       CALL conf_dat2d(dlon_ini, dlat_ini, dlon, dlat, champ)
       CALL inter_barxy(dlon, dlat(:jmdep-1), champ, rlonu(:iim), rlatv, &
            champtime(:, :, l) )
    ENDDO

    call NF95_CLOSE(ncid)

    deallocate(dlon_ini, dlat_ini)
    allocate(yder(lmdep))

    ! interpolation temporelle
    DO j = 1, jjm + 1
       DO i = 1, iim
          yder(:) = SPLINE(timeyear, champtime(i, j, :))
          DO k = 1, 360
             champan(i, j, k) = SPLINT(timeyear, champtime(i, j, :), yder, &
                  real(k-1))
          ENDDO
       ENDDO
    ENDDO
    deallocate(timeyear)

    champan(iim + 1, :, :) = champan(1, :, :)
    forall (k = 1:360) phy_alb(:, k) = pack(champan(:, :, k), dyn_phy)

    DO k = 1, 360
       DO i = 1, klon
          phy_bil(i, k) = 0.0
       ENDDO
    ENDDO

    PRINT *, 'Ecriture du fichier limit'

    call NF95_CREATE("limit.nc", NF90_CLOBBER, nid)

    ierr = NF90_PUT_ATT(nid, NF90_GLOBAL, "title", &
         "Fichier conditions aux limites")
    call NF95_DEF_DIM (nid, "points_physiques", klon, ndim)
    call NF95_DEF_DIM (nid, "time", NF90_UNLIMITED, ntim)

    dims(1) = ndim
    dims(2) = ntim

    ierr = NF90_DEF_VAR (nid, "TEMPS", NF90_FLOAT, ntim, id_tim)
    ierr = NF90_PUT_ATT (nid, id_tim, "title",  &
         "Jour dans l annee")
    ierr = NF90_DEF_VAR (nid, "FOCE", NF90_FLOAT, dims, id_FOCE)
    ierr = NF90_PUT_ATT (nid, id_FOCE, "title", &
         "Fraction ocean")

    ierr = NF90_DEF_VAR (nid, "FSIC", NF90_FLOAT, dims, id_FSIC)
    ierr = NF90_PUT_ATT (nid, id_FSIC, "title", &
         "Fraction glace de mer")

    ierr = NF90_DEF_VAR (nid, "FTER", NF90_FLOAT, dims, id_FTER)
    ierr = NF90_PUT_ATT (nid, id_FTER, "title", &
         "Fraction terre")

    ierr = NF90_DEF_VAR (nid, "FLIC", NF90_FLOAT, dims, id_FLIC)
    ierr = NF90_PUT_ATT (nid, id_FLIC, "title", &
         "Fraction land ice")

    ierr = NF90_DEF_VAR (nid, "SST", NF90_FLOAT, dims, id_SST)
    ierr = NF90_PUT_ATT (nid, id_SST, "title",  &
         "Temperature superficielle de la mer")
    ierr = NF90_DEF_VAR (nid, "BILS", NF90_FLOAT, dims, id_BILS)
    ierr = NF90_PUT_ATT (nid, id_BILS, "title", &
         "Reference flux de chaleur au sol")
    ierr = NF90_DEF_VAR (nid, "ALB", NF90_FLOAT, dims, id_ALB)
    ierr = NF90_PUT_ATT (nid, id_ALB, "title", &
         "Albedo a la surface")
    ierr = NF90_DEF_VAR (nid, "RUG", NF90_FLOAT, dims, id_RUG)
    ierr = NF90_PUT_ATT (nid, id_RUG, "title", &
         "Rugosite")

    call NF95_ENDDEF(nid)

    DO k = 1, 360
       debut(1) = 1
       debut(2) = k

       ierr = NF90_PUT_VAR(nid, id_tim, FLOAT(k), (/k/))
       ierr = NF90_PUT_VAR(nid, id_FOCE, pctsrf_t(:, is_oce, k), debut)
       ierr = NF90_PUT_VAR (nid, id_FSIC, pctsrf_t(:, is_sic, k), debut)
       ierr = NF90_PUT_VAR (nid, id_FTER, pctsrf_t(:, is_ter, k), debut)
       ierr = NF90_PUT_VAR (nid, id_FLIC, pctsrf_t(:, is_lic, k), debut)
       ierr = NF90_PUT_VAR (nid, id_SST, phy_sst(:, k), debut)
       ierr = NF90_PUT_VAR (nid, id_BILS, phy_bil(:, k), debut)
       ierr = NF90_PUT_VAR (nid, id_ALB, phy_alb(:, k), debut)
       ierr = NF90_PUT_VAR (nid, id_RUG, phy_rug(:, k), debut)
    ENDDO

    call NF95_CLOSE(nid)

  END SUBROUTINE limit

end module limit_mod
