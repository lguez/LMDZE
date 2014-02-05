module phyredem_m

  IMPLICIT NONE

contains

  SUBROUTINE phyredem(fichnom, rlat, rlon, pctsrf, tsol, tsoil, tslab, &
       seaice, qsurf, qsol, snow, albedo, alblw, evap, rain_fall, snow_fall, &
       solsw, sollw, fder, radsol, frugs, agesno, zmea, zstd, zsig, zgam, &
       zthe, zpic, zval, t_ancien, q_ancien, rnebcon, ratqs, clwcon, &
       run_off_lic_0, sig1, w01)

    ! From phylmd/phyredem.F, version 1.3 2005/05/25 13:10:09
    ! Author: Z. X. Li (LMD/CNRS)
    ! Date: 1993/08/18
    ! Objet : écriture de l'état de démarrage ou redémarrage pour la physique

    USE dimphy, ONLY: klev, klon, zmasq
    USE dimsoil, ONLY: nsoilmx
    USE indicesol, ONLY: is_lic, is_oce, is_sic, is_ter, nbsrf
    USE netcdf, ONLY: nf90_clobber, nf90_global, nf90_float
    USE netcdf95, ONLY: nf95_create, nf95_put_att, nf95_def_dim, &
         nf95_def_var, nf95_enddef, nf95_redef, nf95_put_var, nf95_close
    USE temps, ONLY: itau_phy

    CHARACTER(len=*) fichnom
    REAL, INTENT(IN):: rlat(klon), rlon(klon)
    REAL, INTENT(IN):: pctsrf(klon, nbsrf)
    REAL tsol(klon, nbsrf)
    REAL tsoil(klon, nsoilmx, nbsrf)
    REAL tslab(klon), seaice(klon) !IM "slab" ocean
    REAL qsurf(klon, nbsrf)
    REAL, intent(in):: qsol(klon)
    REAL snow(klon, nbsrf)
    REAL albedo(klon, nbsrf)
    REAL alblw(klon, nbsrf)
    REAL evap(klon, nbsrf)
    REAL, INTENT(IN):: rain_fall(klon)
    REAL snow_fall(klon)
    REAL solsw(klon)
    REAL, INTENT(IN):: sollw(klon)
    REAL fder(klon)
    REAL radsol(klon)
    REAL frugs(klon, nbsrf)
    REAL agesno(klon, nbsrf)
    REAL, INTENT(IN):: zmea(klon)
    REAL, intent(in):: zstd(klon)
    REAL, intent(in):: zsig(klon)
    REAL zgam(klon)
    REAL zthe(klon)
    REAL zpic(klon)
    REAL zval(klon)
    REAL t_ancien(klon, klev), q_ancien(klon, klev)
    REAL rnebcon(klon, klev), ratqs(klon, klev), clwcon(klon, klev)
    REAL run_off_lic_0(klon)
    real, intent(in):: sig1(klon, klev) ! section adiabatic updraft

    real, intent(in):: w01(klon, klev) 
    ! vertical velocity within adiabatic updraft

    ! Local:

    INTEGER ncid, idim2, idim3
    integer varid, varid_run_off_lic_0, varid_sig1, varid_w01, varid_rlon
    integer varid_rlat, varid_zmasq, varid_fter, varid_flic, varid_foce
    integer varid_fsic
    INTEGER isoil, nsrf
    CHARACTER(len=7) str7
    CHARACTER(len=2) str2

    !------------------------------------------------------------

    PRINT *, 'Call sequence information: phyredem'
    CALL nf95_create(fichnom, nf90_clobber, ncid)

    call nf95_put_att(ncid, nf90_global, 'title', &
         'Fichier redémarrage physique')
    call nf95_put_att(ncid, nf90_global, "itau_phy", itau_phy)

    call nf95_def_dim(ncid, 'points_physiques', klon, idim2)
    call nf95_def_dim(ncid, 'klev', klev, idim3)

    call nf95_def_var(ncid, 'longitude', nf90_float, idim2, varid_rlon)
    call nf95_def_var(ncid, 'latitude', nf90_float, idim2, varid_rlat)

    call nf95_def_var(ncid, 'masque', nf90_float, idim2, varid_zmasq)
    call nf95_put_att(ncid, varid_zmasq, 'title', 'masque terre mer')

    ! Fractions de chaque sous-surface

    call nf95_def_var(ncid, 'FTER', nf90_float, idim2, varid_fter)
    call nf95_put_att(ncid, varid_fter, 'title', 'fraction de continent')

    call nf95_def_var(ncid, 'FLIC', nf90_float, idim2, varid_flic)
    call nf95_put_att(ncid, varid_flic, 'title', 'fraction glace de terre')

    call nf95_def_var(ncid, 'FOCE', nf90_float, idim2, varid_foce)
    call nf95_put_att(ncid, varid_foce, 'title', 'fraction ocean')

    call nf95_def_var(ncid, 'FSIC', nf90_float, idim2, varid_fsic)
    call nf95_put_att(ncid, varid_fsic, 'title', 'fraction glace mer')

    call nf95_enddef(ncid)

    call nf95_put_var(ncid, varid_rlon, rlon)
    call nf95_put_var(ncid, varid_rlat, rlat)
    call nf95_put_var(ncid, varid_zmasq, zmasq)
    call nf95_put_var(ncid, varid_fter, pctsrf(:, is_ter))
    call nf95_put_var(ncid, varid_flic, pctsrf(:, is_lic))
    call nf95_put_var(ncid, varid_foce, pctsrf(:, is_oce))
    call nf95_put_var(ncid, varid_fsic, pctsrf(:, is_sic))

    DO nsrf = 1, nbsrf
       IF (nsrf<=99) THEN
          WRITE (str2, '(i2.2)') nsrf
          call nf95_redef(ncid)
          call nf95_def_var(ncid, 'TS'//str2, nf90_float, idim2, varid)
          call nf95_put_att(ncid, varid, 'title', &
               'Temperature de surface No.'//str2)
          call nf95_enddef(ncid)
       ELSE
          PRINT *, 'Trop de sous-mailles'
          STOP 1
       END IF
       call nf95_put_var(ncid, varid, tsol(:, nsrf))
    END DO

    DO nsrf = 1, nbsrf
       DO isoil = 1, nsoilmx
          IF (isoil<=99 .AND. nsrf<=99) THEN
             WRITE (str7, '(i2.2, "srf", i2.2)') isoil, nsrf
             call nf95_redef(ncid)
             call nf95_def_var(ncid, 'Tsoil'//str7, nf90_float, idim2, varid)
             call nf95_put_att(ncid, varid, 'title', &
                  'Temperature du sol No.'//str7)
             call nf95_enddef(ncid)
          ELSE
             PRINT *, 'Trop de couches'
             STOP 1
          END IF
          call nf95_put_var(ncid, varid, tsoil(:, isoil, nsrf))
       END DO
    END DO

    !IM "slab" ocean
    call nf95_redef(ncid)
    call nf95_def_var(ncid, 'TSLAB', nf90_float, idim2, varid)
    call nf95_put_att(ncid, varid, 'title', &
         'Ecart de la SST (pour slab-ocean)')
    call nf95_enddef(ncid)
    call nf95_put_var(ncid, varid, tslab)

    call nf95_redef(ncid)
    call nf95_def_var(ncid, 'SEAICE', nf90_float, idim2, varid)
    call nf95_put_att(ncid, varid, 'title', &
         'Glace de mer kg/m2 (pour slab-ocean)')
    call nf95_enddef(ncid)
    call nf95_put_var(ncid, varid, seaice)

    DO nsrf = 1, nbsrf
       IF (nsrf<=99) THEN
          WRITE (str2, '(i2.2)') nsrf
          call nf95_redef(ncid)
          call nf95_def_var(ncid, 'QS'//str2, nf90_float, idim2, varid)
          call nf95_put_att(ncid, varid, 'title', &
               'Humidite de surface No.'//str2)
          call nf95_enddef(ncid)
       ELSE
          PRINT *, 'Trop de sous-mailles'
          STOP 1
       END IF
       call nf95_put_var(ncid, varid, qsurf(:, nsrf))
    END DO

    call nf95_redef(ncid)
    call nf95_def_var(ncid, 'QSOL', nf90_float, idim2, varid)
    call nf95_put_att(ncid, varid, 'title', 'Eau dans le sol (mm)')
    call nf95_enddef(ncid)
    call nf95_put_var(ncid, varid, qsol)

    DO nsrf = 1, nbsrf
       IF (nsrf<=99) THEN
          WRITE (str2, '(i2.2)') nsrf
          call nf95_redef(ncid)
          call nf95_def_var(ncid, 'ALBE'//str2, nf90_float, idim2, varid)
          call nf95_put_att(ncid, varid, 'title', &
               'albedo de surface No.'//str2)
          call nf95_enddef(ncid)
       ELSE
          PRINT *, 'Trop de sous-mailles'
          STOP 1
       END IF
       call nf95_put_var(ncid, varid, albedo(:, nsrf))
    END DO

    !IM BEG albedo LW
    DO nsrf = 1, nbsrf
       IF (nsrf<=99) THEN
          WRITE (str2, '(i2.2)') nsrf
          call nf95_redef(ncid)
          call nf95_def_var(ncid, 'ALBLW'//str2, nf90_float, idim2, varid)
          call nf95_put_att(ncid, varid, 'title', &
               'albedo LW de surface No.'//str2)
          call nf95_enddef(ncid)
       ELSE
          PRINT *, 'Trop de sous-mailles'
          STOP 1
       END IF
       call nf95_put_var(ncid, varid, alblw(:, nsrf))
    END DO
    !IM END albedo LW

    DO nsrf = 1, nbsrf
       IF (nsrf<=99) THEN
          WRITE (str2, '(i2.2)') nsrf
          call nf95_redef(ncid)
          call nf95_def_var(ncid, 'EVAP'//str2, nf90_float, idim2, varid)
          call nf95_put_att(ncid, varid, 'title', &
               'Evaporation de surface No.'//str2)
          call nf95_enddef(ncid)
       ELSE
          PRINT *, 'Trop de sous-mailles'
          STOP 1
       END IF
       call nf95_put_var(ncid, varid, evap(:, nsrf))
    END DO

    DO nsrf = 1, nbsrf
       IF (nsrf<=99) THEN
          WRITE (str2, '(i2.2)') nsrf
          call nf95_redef(ncid)
          call nf95_def_var(ncid, 'SNOW'//str2, nf90_float, idim2, varid)
          call nf95_put_att(ncid, varid, 'title', &
               'Neige de surface No.'//str2)
          call nf95_enddef(ncid)
       ELSE
          PRINT *, 'Trop de sous-mailles'
          STOP 1
       END IF
       call nf95_put_var(ncid, varid, snow(:, nsrf))
    END DO

    call nf95_redef(ncid)
    call nf95_def_var(ncid, 'RADS', nf90_float, idim2, varid)
    call nf95_put_att(ncid, varid, 'title', &
         'Rayonnement net a la surface')
    call nf95_enddef(ncid)
    call nf95_put_var(ncid, varid, radsol)

    call nf95_redef(ncid)
    call nf95_def_var(ncid, 'solsw', nf90_float, idim2, varid)
    call nf95_put_att(ncid, varid, 'title', &
         'Rayonnement solaire a la surface')
    call nf95_enddef(ncid)
    call nf95_put_var(ncid, varid, solsw)

    call nf95_redef(ncid)
    call nf95_def_var(ncid, 'sollw', nf90_float, idim2, varid)
    call nf95_put_att(ncid, varid, 'title', &
         'Rayonnement IF a la surface')
    call nf95_enddef(ncid)
    call nf95_put_var(ncid, varid, sollw)

    call nf95_redef(ncid)
    call nf95_def_var(ncid, 'fder', nf90_float, idim2, varid)
    call nf95_put_att(ncid, varid, 'title', 'Derive de flux')
    call nf95_enddef(ncid)
    call nf95_put_var(ncid, varid, fder)

    call nf95_redef(ncid)
    call nf95_def_var(ncid, 'rain_f', nf90_float, idim2, varid)
    call nf95_put_att(ncid, varid, 'title', 'precipitation liquide')
    call nf95_enddef(ncid)
    call nf95_put_var(ncid, varid, rain_fall)

    call nf95_redef(ncid)
    call nf95_def_var(ncid, 'snow_f', nf90_float, idim2, varid)
    call nf95_put_att(ncid, varid, 'title', 'precipitation solide')
    call nf95_enddef(ncid)
    call nf95_put_var(ncid, varid, snow_fall)

    DO nsrf = 1, nbsrf
       IF (nsrf<=99) THEN
          WRITE (str2, '(i2.2)') nsrf
          call nf95_redef(ncid)
          call nf95_def_var(ncid, 'RUG'//str2, nf90_float, idim2, varid)
          call nf95_put_att(ncid, varid, 'title', &
               'rugosite de surface No.'//str2)
          call nf95_enddef(ncid)
       ELSE
          PRINT *, 'Trop de sous-mailles'
          STOP 1
       END IF
       call nf95_put_var(ncid, varid, frugs(:, nsrf))
    END DO

    DO nsrf = 1, nbsrf
       IF (nsrf<=99) THEN
          WRITE (str2, '(i2.2)') nsrf
          call nf95_redef(ncid)
          call nf95_def_var(ncid, 'AGESNO'//str2, nf90_float, idim2, varid)
          call nf95_put_att(ncid, varid, 'title', &
               'Age de la neige surface No.'//str2)
          call nf95_enddef(ncid)
       ELSE
          PRINT *, 'Trop de sous-mailles'
          STOP 1
       END IF
       call nf95_put_var(ncid, varid, agesno(:, nsrf))
    END DO

    call nf95_redef(ncid)
    call nf95_def_var(ncid, 'ZMEA', nf90_float, idim2, varid)
    call nf95_enddef(ncid)
    call nf95_put_var(ncid, varid, zmea)

    call nf95_redef(ncid)
    call nf95_def_var(ncid, 'ZSTD', nf90_float, idim2, varid)
    call nf95_enddef(ncid)
    call nf95_put_var(ncid, varid, zstd)
    call nf95_redef(ncid)
    call nf95_def_var(ncid, 'ZSIG', nf90_float, idim2, varid)
    call nf95_enddef(ncid)
    call nf95_put_var(ncid, varid, zsig)
    call nf95_redef(ncid)
    call nf95_def_var(ncid, 'ZGAM', nf90_float, idim2, varid)
    call nf95_enddef(ncid)
    call nf95_put_var(ncid, varid, zgam)
    call nf95_redef(ncid)
    call nf95_def_var(ncid, 'ZTHE', nf90_float, idim2, varid)
    call nf95_enddef(ncid)
    call nf95_put_var(ncid, varid, zthe)
    call nf95_redef(ncid)
    call nf95_def_var(ncid, 'ZPIC', nf90_float, idim2, varid)
    call nf95_enddef(ncid)
    call nf95_put_var(ncid, varid, zpic)
    call nf95_redef(ncid)
    call nf95_def_var(ncid, 'ZVAL', nf90_float, idim2, varid)
    call nf95_enddef(ncid)
    call nf95_put_var(ncid, varid, zval)

    call nf95_redef(ncid)
    call nf95_def_var(ncid, 'TANCIEN', nf90_float, (/idim2, idim3/), varid)
    call nf95_enddef(ncid)
    call nf95_put_var(ncid, varid, t_ancien)

    call nf95_redef(ncid)
    call nf95_def_var(ncid, 'QANCIEN', nf90_float, (/idim2, idim3/), varid)
    call nf95_enddef(ncid)
    call nf95_put_var(ncid, varid, q_ancien)

    call nf95_redef(ncid)
    call nf95_def_var(ncid, 'RUGMER', nf90_float, idim2, varid)
    call nf95_put_att(ncid, varid, 'title', &
         'Longueur de rugosite sur mer')
    call nf95_enddef(ncid)
    call nf95_put_var(ncid, varid, frugs(:, is_oce))

    call nf95_redef(ncid)
    call nf95_def_var(ncid, 'CLWCON', nf90_float, idim2, varid)
    call nf95_put_att(ncid, varid, 'title', 'Eau liquide convective')
    call nf95_enddef(ncid)
    call nf95_put_var(ncid, varid, clwcon(:, 1))

    call nf95_redef(ncid)
    call nf95_def_var(ncid, 'RNEBCON', nf90_float, idim2, varid)
    call nf95_put_att(ncid, varid, 'title', 'Nebulosite convective')
    call nf95_enddef(ncid)
    call nf95_put_var(ncid, varid, rnebcon(:, 1))

    call nf95_redef(ncid)
    call nf95_def_var(ncid, 'RATQS', nf90_float, idim2, varid)
    call nf95_put_att(ncid, varid, 'title', 'Ratqs')
    call nf95_enddef(ncid)
    call nf95_put_var(ncid, varid, ratqs(:, 1))

    call nf95_redef(ncid)
    call nf95_def_var(ncid, 'RUNOFFLIC0', nf90_float, idim2, &
         varid_run_off_lic_0)
    call nf95_put_att(ncid, varid_run_off_lic_0, 'title', 'Runofflic0')

    call nf95_def_var(ncid, 'sig1', nf90_float, (/idim2, idim3/), varid_sig1)
    call nf95_put_att(ncid, varid_sig1, 'long_name', &
         'section adiabatic updraft')

    call nf95_def_var(ncid, 'w01', nf90_float, (/idim2, idim3/), varid_w01)
    call nf95_put_att(ncid, varid_w01, 'long_name', &
         'vertical velocity within adiabatic updraft')

    call nf95_enddef(ncid)

    call nf95_put_var(ncid, varid_run_off_lic_0, run_off_lic_0)
    call nf95_put_var(ncid, varid_sig1, sig1)
    call nf95_put_var(ncid, varid_w01, w01)

    call nf95_close(ncid)

  END SUBROUTINE phyredem

end module phyredem_m
