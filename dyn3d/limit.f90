module limit_mod

  IMPLICIT none

contains

  SUBROUTINE limit(pctsrf)

    ! Authors: L. Fairhead, Z. X. Li, P. Le Van

    ! This subroutine creates files containing boundary conditions.
    ! It uses files with climatological data.  Both grids must be
    ! regular.

    use conf_dat2d_m, only: conf_dat2d
    use dimensions, only: iim, jjm
    use dimphy, only: klon
    use dynetat0_m, only: rlonu, rlatv
    use grid_change, only: dyn_phy
    use indicesol, only: epsfra, is_ter, is_oce, is_lic, is_sic
    use inter_barxy_m, only: inter_barxy
    use netcdf95, only: NF95_CLOSE, NF95_CREATE, NF95_DEF_DIM, nf95_def_var, &
         nf95_enddef, nf95_get_var, nf95_gw_var, nf95_inq_dimid, &
         nf95_inq_varid, nf95_inquire_dimension, NF95_OPEN, NF95_PUT_ATT, &
         NF95_PUT_VAR
    use netcdf, only: NF90_CLOBBER, NF90_FLOAT, NF90_GLOBAL, NF90_NOWRITE, &
         NF90_UNLIMITED
    use nr_util, only: assert
    use numer_rec_95, only: spline, splint
    use phyetat0_m, only: masque
    use start_init_orog_m, only: mask
    use unit_nml_m, only: unit_nml

    REAL, intent(inout):: pctsrf(:, :) ! (klon, nbsrf)
    ! "pctsrf(i, :)" is the composition of the surface at horizontal
    ! position "i".

    ! Local:

    LOGICAL:: extrap = .FALSE.
    ! (extrapolation de donn\'ees, comme pour les SST lorsque le fichier
    ! ne contient pas uniquement des points oc\'eaniques)

    REAL phy_bil(klon, 360)
    REAL phy_ice(klon)

    ! Pour le champ de d\'epart:
    INTEGER imdep, jmdep, lmdep

    REAL, ALLOCATABLE:: dlon(:), dlat(:)
    REAL, ALLOCATABLE:: dlon_ini(:), dlat_ini(:), timeyear(:)
    REAL, ALLOCATABLE:: champ(:, :)
    REAL, ALLOCATABLE:: work(:, :)

    ! Pour le champ interpol\'e 3D :
    REAL, allocatable:: champtime(:, :, :)
    REAL champan(iim + 1, jjm + 1, 360)

    ! Pour l'inteprolation verticale :
    REAL, allocatable:: yder(:)

    INTEGER ndim, ntim
    INTEGER varid_time
    INTEGER id_SST, id_BILS, id_RUG, id_ALB
    INTEGER id_FOCE, id_FSIC, id_FTER, id_FLIC

    INTEGER i, j, k, l
    INTEGER ncid, ncid_limit, varid, dimid

    REAL, parameter:: tmidmonth(12) = [(15. + 30. * i, i = 0, 11)]

    namelist /limit_nml/extrap

    !--------------------

    print *, "Call sequence information: limit"

    print *, "Enter namelist 'limit_nml'."
    read(unit=*, nml=limit_nml)
    write(unit_nml, nml=limit_nml)

    call NF95_CREATE("limit.nc", NF90_CLOBBER, ncid_limit)

    call NF95_PUT_ATT(ncid_limit, NF90_GLOBAL, "title", &
         "Fichier conditions aux limites")
    call NF95_DEF_DIM(ncid_limit, "points_physiques", klon, ndim)
    call NF95_DEF_DIM(ncid_limit, "time", NF90_UNLIMITED, ntim)

    call NF95_DEF_VAR(ncid_limit, "TEMPS", NF90_FLOAT, ntim, varid_time)
    call NF95_PUT_ATT(ncid_limit, varid_time, "title", "Jour dans l annee")

    call NF95_DEF_VAR(ncid_limit, "FOCE", NF90_FLOAT, dimids=[ndim, ntim], &
         varid=id_foce)
    call NF95_PUT_ATT(ncid_limit, id_FOCE, "title", "Fraction ocean")

    call NF95_DEF_VAR(ncid_limit, "FSIC", NF90_FLOAT, dimids=[ndim, ntim], &
         varid=id_FSIC)
    call NF95_PUT_ATT(ncid_limit, id_FSIC, "title", "Fraction glace de mer")

    call NF95_DEF_VAR(ncid_limit, "FTER", NF90_FLOAT, dimids=[ndim, ntim], &
         varid=id_FTER)
    call NF95_PUT_ATT(ncid_limit, id_FTER, "title", "Fraction terre")

    call NF95_DEF_VAR(ncid_limit, "FLIC", NF90_FLOAT, dimids=[ndim, ntim], &
         varid=id_FLIC)
    call NF95_PUT_ATT(ncid_limit, id_FLIC, "title", "Fraction land ice")

    call NF95_DEF_VAR(ncid_limit, "SST", NF90_FLOAT, dimids=[ndim, ntim], &
         varid=id_SST)
    call NF95_PUT_ATT(ncid_limit, id_SST, "title",  &
         "Temperature superficielle de la mer")

    call NF95_DEF_VAR(ncid_limit, "BILS", NF90_FLOAT, dimids=[ndim, ntim], &
         varid=id_BILS)
    call NF95_PUT_ATT(ncid_limit, id_BILS, "title", &
         "Reference flux de chaleur au sol")

    call NF95_DEF_VAR(ncid_limit, "ALB", NF90_FLOAT, dimids=[ndim, ntim], &
         varid=id_ALB)
    call NF95_PUT_ATT(ncid_limit, id_ALB, "title", "Albedo a la surface")

    call NF95_DEF_VAR(ncid_limit, "RUG", NF90_FLOAT, dimids=[ndim, ntim], &
         varid=id_RUG)
    call NF95_PUT_ATT(ncid_limit, id_RUG, "title", "Rugosite")

    call NF95_ENDDEF(ncid_limit)

    call NF95_PUT_VAR(ncid_limit, varid_time, [(k, k = 1, 360)])
    
    PRINT *, 'Processing rugosity...'

    call NF95_OPEN('Rugos.nc', NF90_NOWRITE, ncid)

    ! Read coordinate variables:

    call nf95_inq_varid(ncid, "longitude", varid)
    call nf95_gw_var(ncid, varid, dlon_ini)
    imdep = size(dlon_ini)

    call nf95_inq_varid(ncid, "latitude", varid)
    call nf95_gw_var(ncid, varid, dlat_ini)
    jmdep = size(dlat_ini)

    call nf95_inq_varid(ncid, "temps", varid)
    call nf95_gw_var(ncid, varid, timeyear)
    lmdep = size(timeyear)

    ALLOCATE(champ(imdep, jmdep), champtime(iim, jjm + 1, lmdep))
    allocate(dlon(imdep), dlat(jmdep))
    call NF95_INQ_VARID(ncid, 'RUGOS', varid)

    ! Read the primary variable day by day and regrid horizontally,
    ! result in "champtime":
    DO  l = 1, lmdep
       call NF95_GET_VAR(ncid, varid, champ, start=[1, 1, l])
       CALL conf_dat2d(dlon_ini, dlat_ini, dlon, dlat, champ)
       CALL inter_barxy(dlon, dlat(:jmdep -1), LOG(champ), rlonu(:iim), &
            rlatv, champtime(:, :, l))
       champtime(:, :, l) = EXP(champtime(:, :, l))
       where (nint(mask(:iim, :)) /= 1) champtime(:, :, l) = 0.001
    end do

    call NF95_CLOSE(ncid)

    DEALLOCATE(dlon, dlat, champ)
    allocate(yder(lmdep))

    ! Interpolate monthly values to daily values, at each horizontal position:
    DO j = 1, jjm + 1
       DO i = 1, iim
          yder = SPLINE(timeyear, champtime(i, j, :))
          DO k = 1, 360
             champan(i, j, k) = SPLINT(timeyear, champtime(i, j, :), yder, &
                  real(k-1))
          ENDDO
       ENDDO
    ENDDO

    deallocate(champtime, yder)
    champan(iim + 1, :, :) = champan(1, :, :)

    DO k = 1, 360
       call NF95_PUT_VAR(ncid_limit, id_RUG, pack(champan(:, :, k), dyn_phy), &
            start=[1, k])
    ENDDO

    PRINT *, 'Processing sea ice...'
    call NF95_OPEN('amipbc_sic_1x1.nc', NF90_NOWRITE, ncid)

    call nf95_inq_varid(ncid, "longitude", varid)
    call nf95_gw_var(ncid, varid, dlon_ini)
    imdep = size(dlon_ini)

    call nf95_inq_varid(ncid, "latitude", varid)
    call nf95_gw_var(ncid, varid, dlat_ini)
    jmdep = size(dlat_ini)

    call nf95_inq_dimid(ncid, "time", dimid)
    call NF95_INQuire_DIMension(ncid, dimid, nclen=lmdep)
    print *, 'lmdep = ', lmdep
    ! Coordonn\'ee temporelle fichiers AMIP pas en jours. Ici on suppose
    ! qu'on a 12 mois (de 30 jours).
    IF (lmdep /= 12) then
       print *, 'Unknown AMIP file: not 12 months?'
       STOP 1
    end IF

    ALLOCATE(champ(imdep, jmdep), champtime(iim, jjm + 1, lmdep))
    ALLOCATE(dlon(imdep), dlat(jmdep))
    call NF95_INQ_VARID(ncid, 'sicbcs', varid)
    DO l = 1, lmdep
       call NF95_GET_VAR(ncid, varid, champ, start=[1, 1, l])
       CALL conf_dat2d(dlon_ini, dlat_ini, dlon, dlat, champ)
       CALL inter_barxy(dlon, dlat(:jmdep -1), champ, rlonu(:iim), rlatv, &
            champtime(:, :, l))
    ENDDO

    call NF95_CLOSE(ncid)

    DEALLOCATE(dlon, dlat, champ)
    PRINT *, 'Interpolation temporelle'
    allocate(yder(lmdep))

    DO j = 1, jjm + 1
       DO i = 1, iim
          yder = SPLINE(tmidmonth, champtime(i, j, :))
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
       phy_ice = pack(champan(:, :, k), dyn_phy)

       ! (utilisation de la sous-maille fractionnelle tandis que l'ancien
       ! codage utilisait l'indicateur du sol (0, 1, 2, 3))
       ! PB en attendant de mettre fraction de terre
       WHERE (phy_ice < EPSFRA) phy_ice = 0.

       pctsrf(:, is_sic) = max(phy_ice - pctsrf(:, is_lic), 0.)
       ! Il y a des cas o\`u il y a de la glace dans landiceref et
       ! pas dans AMIP
       WHERE (1. - masque < EPSFRA)
          pctsrf(:, is_sic) = 0.
          pctsrf(:, is_oce) = 0.
       elsewhere
          where (pctsrf(:, is_sic) >= 1 - masque)
             pctsrf(:, is_sic) = 1. - masque
             pctsrf(:, is_oce) = 0.
          ELSEwhere
             pctsrf(:, is_oce) = 1. - masque - pctsrf(:, is_sic)
             where (pctsrf(:, is_oce) < EPSFRA)
                pctsrf(:, is_oce) = 0.
                pctsrf(:, is_sic) = 1 - masque
             end where
          end where
       end where

       DO i = 1, klon
          if (pctsrf(i, is_oce) < 0.) then
             print *, "k = ", k
             print *, 'Bad surface fraction: pctsrf(', i, ', is_oce) = ', &
                  pctsrf(i, is_oce)
          ENDIF
          IF (abs(sum(pctsrf(i, :)) - 1.) > EPSFRA) THEN 
             print *, "k = ", k
             print *, 'Bad surface fraction:'
             print *, "pctsrf(", i, ", :) = ", pctsrf(i, :)
             print *, "phy_ice(", i, ") = ", phy_ice(i)
          ENDIF
       END DO

       call NF95_PUT_VAR(ncid_limit, id_FOCE, pctsrf(:, is_oce), start=[1, k])
       call NF95_PUT_VAR(ncid_limit, id_FSIC, pctsrf(:, is_sic), start=[1, k])
       call NF95_PUT_VAR(ncid_limit, id_FTER, pctsrf(:, is_ter), start=[1, k])
       call NF95_PUT_VAR(ncid_limit, id_FLIC, pctsrf(:, is_lic), start=[1, k])
    end DO
    
    PRINT *, 'Traitement de la sst'
    call NF95_OPEN('amipbc_sst_1x1.nc', NF90_NOWRITE, ncid)

    call nf95_inq_varid(ncid, "longitude", varid)
    call nf95_gw_var(ncid, varid, dlon_ini)
    imdep = size(dlon_ini)

    call nf95_inq_varid(ncid, "latitude", varid)
    call nf95_gw_var(ncid, varid, dlat_ini)
    jmdep = size(dlat_ini)

    call nf95_inq_dimid(ncid, "time", dimid)
    call NF95_INQuire_DIMension(ncid, dimid, nclen=lmdep)
    print *, 'lmdep = ', lmdep
    ! Ici on suppose qu'on a 12 mois (de 30 jours).
    call assert(lmdep == 12, 'limit: AMIP file does not contain 12 months')

    ALLOCATE(champ(imdep, jmdep), champtime(iim, jjm + 1, lmdep))
    IF (extrap) ALLOCATE(work(imdep, jmdep))
    ALLOCATE(dlon(imdep), dlat(jmdep))
    call NF95_INQ_VARID(ncid, 'tosbcs', varid)

    DO l = 1, lmdep
       call NF95_GET_VAR(ncid, varid, champ, start=[1, 1, l])
       CALL conf_dat2d(dlon_ini, dlat_ini, dlon, dlat, champ)
       IF (extrap) &
            CALL extrapol(champ, imdep, jmdep, 999999., .TRUE., .TRUE., 2, work)
       CALL inter_barxy(dlon, dlat(:jmdep -1), champ, rlonu(:iim), rlatv, &
            champtime(:, :, l))
    ENDDO

    call NF95_CLOSE(ncid)
    DEALLOCATE(dlon, dlat, champ)
    allocate(yder(lmdep))

    ! interpolation temporelle
    DO j = 1, jjm + 1
       DO i = 1, iim
          yder = SPLINE(tmidmonth, champtime(i, j, :))
          DO k = 1, 360
             champan(i, j, k) = SPLINT(tmidmonth, champtime(i, j, :), yder, &
                  real(k-1))
          ENDDO
       ENDDO
    ENDDO

    deallocate(champtime, yder)
    champan(iim + 1, :, :) = champan(1, :, :)

    !IM14/03/2002 : SST amipbc greater then 271.38
    PRINT *, 'limit: SST Amipbc >= 271.38 '

    DO k = 1, 360
       DO j = 1, jjm + 1
          DO i = 1, iim
             champan(i, j, k) = max(champan(i, j, k), 271.38)
          ENDDO
          
          champan(iim + 1, j, k) = champan(1, j, k)
       ENDDO
    ENDDO
    
    DO k = 1, 360
       call NF95_PUT_VAR(ncid_limit, id_SST, pack(champan(:, :, k), dyn_phy), &
            start=[1, k])
    end DO

    PRINT *, "Traitement de l'albedo..."
    call NF95_OPEN('Albedo.nc', NF90_NOWRITE, ncid)

    call nf95_inq_varid(ncid, "longitude", varid)
    call nf95_gw_var(ncid, varid, dlon_ini)
    imdep = size(dlon_ini)

    call nf95_inq_varid(ncid, "latitude", varid)
    call nf95_gw_var(ncid, varid, dlat_ini)
    jmdep = size(dlat_ini)

    call nf95_inq_varid(ncid, "temps", varid)
    call nf95_gw_var(ncid, varid, timeyear)
    lmdep = size(timeyear)

    ALLOCATE(champ(imdep, jmdep), champtime(iim, jjm + 1, lmdep))
    ALLOCATE(dlon(imdep), dlat(jmdep))
    call NF95_INQ_VARID(ncid, 'ALBEDO', varid)

    DO l = 1, lmdep
       PRINT *, "timeyear(", l, ") =", timeyear(l)
       call NF95_GET_VAR(ncid, varid, champ, start=[1, 1, l])
       CALL conf_dat2d(dlon_ini, dlat_ini, dlon, dlat, champ)
       CALL inter_barxy(dlon, dlat(:jmdep-1), champ, rlonu(:iim), rlatv, &
            champtime(:, :, l))
    ENDDO

    call NF95_CLOSE(ncid)

    allocate(yder(lmdep))

    ! interpolation temporelle
    DO j = 1, jjm + 1
       DO i = 1, iim
          yder = SPLINE(timeyear, champtime(i, j, :))
          DO k = 1, 360
             champan(i, j, k) = SPLINT(timeyear, champtime(i, j, :), yder, &
                  real(k-1))
          ENDDO
       ENDDO
    ENDDO

    champan(iim + 1, :, :) = champan(1, :, :)

    DO k = 1, 360
       call NF95_PUT_VAR(ncid_limit, id_ALB, pack(champan(:, :, k), dyn_phy), &
            start=[1, k])
    end DO

    phy_bil = 0.
    call NF95_PUT_VAR(ncid_limit, id_BILS, phy_bil)

    call NF95_CLOSE(ncid_limit)

  END SUBROUTINE limit

end module limit_mod
