module phyredem_m

  IMPLICIT NONE

contains

  SUBROUTINE phyredem(fichnom, rlat, rlon, pctsrf, tsol, tsoil, tslab, &
       seaice, qsurf, qsol, snow, albedo, alblw, evap, rain_fall,&
       snow_fall, solsw, sollw, fder, radsol, frugs, agesno, zmea,&
       zstd, zsig, zgam, zthe, zpic, zval, t_ancien, q_ancien,&
       rnebcon, ratqs, clwcon, run_off_lic_0)

    ! From phylmd/phyredem.F, v 1.3 2005/05/25 13:10:09
    ! Auteur(s) Z.X. Li (LMD/CNRS) date: 19930818
    ! Objet: Ecriture de l'etat de démarrage ou redémarrage pour la physique

    USE indicesol, ONLY : is_lic, is_oce, is_sic, is_ter, nbsrf
    USE dimphy, ONLY : klev, klon, zmasq
    USE dimsoil, ONLY : nsoilmx
    USE temps, ONLY : itau_phy
    USE netcdf95, ONLY : nf95_create, nf95_put_att, nf95_def_dim, &
         nf95_def_var, nf95_enddef, nf95_redef, nf95_put_var, nf95_close
    USE netcdf, ONLY : nf90_clobber, nf90_global, nf90_float

    CHARACTER(len=*) fichnom
    REAL, INTENT (IN) :: rlat(klon), rlon(klon)
    REAL :: tsol(klon, nbsrf)
    REAL :: tsoil(klon, nsoilmx, nbsrf)

    REAL :: tslab(klon), seaice(klon) !IM "slab" ocean
    REAL :: qsurf(klon, nbsrf)
    REAL :: qsol(klon)
    REAL :: snow(klon, nbsrf)
    REAL :: albedo(klon, nbsrf)

    REAL :: alblw(klon, nbsrf)

    REAL :: evap(klon, nbsrf)
    REAL :: rain_fall(klon)
    REAL :: snow_fall(klon)
    REAL :: solsw(klon)
    REAL :: sollw(klon)
    REAL :: fder(klon)
    REAL :: radsol(klon)
    REAL :: frugs(klon, nbsrf)
    REAL :: agesno(klon, nbsrf)
    REAL :: zmea(klon)
    REAL, intent(in):: zstd(klon)
    REAL, intent(in):: zsig(klon)
    REAL :: zgam(klon)
    REAL :: zthe(klon)
    REAL :: zpic(klon)
    REAL :: zval(klon)
    REAL :: pctsrf(klon, nbsrf)
    REAL :: t_ancien(klon, klev), q_ancien(klon, klev)
    REAL :: clwcon(klon, klev), rnebcon(klon, klev), ratqs(klon, klev)
    REAL :: run_off_lic_0(klon)

    INTEGER :: nid, nvarid, idim2, idim3

    INTEGER :: isoil, nsrf
    CHARACTER (7) :: str7
    CHARACTER (2) :: str2

    !------------------------------------------------------------

    PRINT *, 'Call sequence information: phyredem'
    CALL nf95_create(fichnom, nf90_clobber, nid)

    call nf95_put_att(nid, nf90_global, 'title', &
         'Fichier redémarrage physique')
    call nf95_put_att(nid, nf90_global, "itau_phy", itau_phy)

    call nf95_def_dim(nid, 'points_physiques', klon, idim2)
    call nf95_def_dim(nid, 'horizon_vertical', klon*klev, idim3)

    call nf95_def_var(nid, 'longitude', nf90_float, idim2, nvarid)
    call nf95_put_att(nid, nvarid, 'title', &
         'Longitudes de la grille physique')
    call nf95_enddef(nid)
    call nf95_put_var(nid, nvarid, rlon)

    call nf95_redef(nid)
    call nf95_def_var(nid, 'latitude', nf90_float, idim2, nvarid)
    call nf95_put_att(nid, nvarid, 'title', &
         'Latitudes de la grille physique')
    call nf95_enddef(nid)
    call nf95_put_var(nid, nvarid, rlat)

    ! PB ajout du masque terre/mer

    call nf95_redef(nid)
    call nf95_def_var(nid, 'masque', nf90_float, idim2, nvarid)
    call nf95_put_att(nid, nvarid, 'title', 'masque terre mer')
    call nf95_enddef(nid)
    call nf95_put_var(nid, nvarid, zmasq)
    ! BP ajout des fraction de chaque sous-surface

    ! 1. fraction de terre

    call nf95_redef(nid)
    call nf95_def_var(nid, 'FTER', nf90_float, idim2, nvarid)
    call nf95_put_att(nid, nvarid, 'title', 'fraction de continent')
    call nf95_enddef(nid)
    call nf95_put_var(nid, nvarid, pctsrf(:, is_ter))

    ! 2. Fraction de glace de terre

    call nf95_redef(nid)
    call nf95_def_var(nid, 'FLIC', nf90_float, idim2, nvarid)
    call nf95_put_att(nid, nvarid, 'title', 'fraction glace de terre')
    call nf95_enddef(nid)
    call nf95_put_var(nid, nvarid, pctsrf(:, is_lic))

    ! 3. fraction ocean

    call nf95_redef(nid)
    call nf95_def_var(nid, 'FOCE', nf90_float, idim2, nvarid)
    call nf95_put_att(nid, nvarid, 'title', 'fraction ocean')
    call nf95_enddef(nid)
    call nf95_put_var(nid, nvarid, pctsrf(:, is_oce))

    ! 4. Fraction glace de mer

    call nf95_redef(nid)
    call nf95_def_var(nid, 'FSIC', nf90_float, idim2, nvarid)
    call nf95_put_att(nid, nvarid, 'title', 'fraction glace mer')
    call nf95_enddef(nid)
    call nf95_put_var(nid, nvarid, pctsrf(:, is_sic))



    DO nsrf = 1, nbsrf
       IF (nsrf<=99) THEN
          WRITE (str2, '(i2.2)') nsrf
          call nf95_redef(nid)
          call nf95_def_var(nid, 'TS'//str2, nf90_float, idim2, nvarid)
          call nf95_put_att(nid, nvarid, 'title', &
               'Temperature de surface No.'//str2)
          call nf95_enddef(nid)
       ELSE
          PRINT *, 'Trop de sous-mailles'
          STOP 1
       END IF
       call nf95_put_var(nid, nvarid, tsol(:, nsrf))
    END DO

    DO nsrf = 1, nbsrf
       DO isoil = 1, nsoilmx
          IF (isoil<=99 .AND. nsrf<=99) THEN
             WRITE (str7, '(i2.2, "srf", i2.2)') isoil, nsrf
             call nf95_redef(nid)
             call nf95_def_var(nid, 'Tsoil'//str7, nf90_float, idim2, nvarid)
             call nf95_put_att(nid, nvarid, 'title', &
                  'Temperature du sol No.'//str7)
             call nf95_enddef(nid)
          ELSE
             PRINT *, 'Trop de couches'
             STOP 1
          END IF
          call nf95_put_var(nid, nvarid, tsoil(:, isoil, nsrf))
       END DO
    END DO

    !IM "slab" ocean
    call nf95_redef(nid)
    call nf95_def_var(nid, 'TSLAB', nf90_float, idim2, nvarid)
    call nf95_put_att(nid, nvarid, 'title', &
         'Ecart de la SST (pour slab-ocean)')
    call nf95_enddef(nid)
    call nf95_put_var(nid, nvarid, tslab)

    call nf95_redef(nid)
    call nf95_def_var(nid, 'SEAICE', nf90_float, idim2, nvarid)
    call nf95_put_att(nid, nvarid, 'title', &
         'Glace de mer kg/m2 (pour slab-ocean)')
    call nf95_enddef(nid)
    call nf95_put_var(nid, nvarid, seaice)

    DO nsrf = 1, nbsrf
       IF (nsrf<=99) THEN
          WRITE (str2, '(i2.2)') nsrf
          call nf95_redef(nid)
          call nf95_def_var(nid, 'QS'//str2, nf90_float, idim2, nvarid)
          call nf95_put_att(nid, nvarid, 'title', &
               'Humidite de surface No.'//str2)
          call nf95_enddef(nid)
       ELSE
          PRINT *, 'Trop de sous-mailles'
          STOP 1
       END IF
       call nf95_put_var(nid, nvarid, qsurf(:, nsrf))
    END DO

    call nf95_redef(nid)
    call nf95_def_var(nid, 'QSOL', nf90_float, idim2, nvarid)
    call nf95_put_att(nid, nvarid, 'title', 'Eau dans le sol (mm)')
    call nf95_enddef(nid)
    call nf95_put_var(nid, nvarid, qsol)

    DO nsrf = 1, nbsrf
       IF (nsrf<=99) THEN
          WRITE (str2, '(i2.2)') nsrf
          call nf95_redef(nid)
          call nf95_def_var(nid, 'ALBE'//str2, nf90_float, idim2, nvarid)
          call nf95_put_att(nid, nvarid, 'title', &
               'albedo de surface No.'//str2)
          call nf95_enddef(nid)
       ELSE
          PRINT *, 'Trop de sous-mailles'
          STOP 1
       END IF
       call nf95_put_var(nid, nvarid, albedo(:, nsrf))
    END DO

    !IM BEG albedo LW
    DO nsrf = 1, nbsrf
       IF (nsrf<=99) THEN
          WRITE (str2, '(i2.2)') nsrf
          call nf95_redef(nid)
          call nf95_def_var(nid, 'ALBLW'//str2, nf90_float, idim2, nvarid)
          call nf95_put_att(nid, nvarid, 'title', &
               'albedo LW de surface No.'//str2)
          call nf95_enddef(nid)
       ELSE
          PRINT *, 'Trop de sous-mailles'
          STOP 1
       END IF
       call nf95_put_var(nid, nvarid, alblw(:, nsrf))
    END DO
    !IM END albedo LW

    DO nsrf = 1, nbsrf
       IF (nsrf<=99) THEN
          WRITE (str2, '(i2.2)') nsrf
          call nf95_redef(nid)
          call nf95_def_var(nid, 'EVAP'//str2, nf90_float, idim2, nvarid)
          call nf95_put_att(nid, nvarid, 'title', &
               'Evaporation de surface No.'//str2)
          call nf95_enddef(nid)
       ELSE
          PRINT *, 'Trop de sous-mailles'
          STOP 1
       END IF
       call nf95_put_var(nid, nvarid, evap(:, nsrf))
    END DO


    DO nsrf = 1, nbsrf
       IF (nsrf<=99) THEN
          WRITE (str2, '(i2.2)') nsrf
          call nf95_redef(nid)
          call nf95_def_var(nid, 'SNOW'//str2, nf90_float, idim2, nvarid)
          call nf95_put_att(nid, nvarid, 'title', &
               'Neige de surface No.'//str2)
          call nf95_enddef(nid)
       ELSE
          PRINT *, 'Trop de sous-mailles'
          STOP 1
       END IF
       call nf95_put_var(nid, nvarid, snow(:, nsrf))
    END DO


    call nf95_redef(nid)
    call nf95_def_var(nid, 'RADS', nf90_float, idim2, nvarid)
    call nf95_put_att(nid, nvarid, 'title', &
         'Rayonnement net a la surface')
    call nf95_enddef(nid)
    call nf95_put_var(nid, nvarid, radsol)

    call nf95_redef(nid)
    call nf95_def_var(nid, 'solsw', nf90_float, idim2, nvarid)
    call nf95_put_att(nid, nvarid, 'title', &
         'Rayonnement solaire a la surface')
    call nf95_enddef(nid)
    call nf95_put_var(nid, nvarid, solsw)

    call nf95_redef(nid)
    call nf95_def_var(nid, 'sollw', nf90_float, idim2, nvarid)
    call nf95_put_att(nid, nvarid, 'title', &
         'Rayonnement IF a la surface')
    call nf95_enddef(nid)
    call nf95_put_var(nid, nvarid, sollw)

    call nf95_redef(nid)
    call nf95_def_var(nid, 'fder', nf90_float, idim2, nvarid)
    call nf95_put_att(nid, nvarid, 'title', 'Derive de flux')
    call nf95_enddef(nid)
    call nf95_put_var(nid, nvarid, fder)

    call nf95_redef(nid)
    call nf95_def_var(nid, 'rain_f', nf90_float, idim2, nvarid)
    call nf95_put_att(nid, nvarid, 'title', 'precipitation liquide')
    call nf95_enddef(nid)
    call nf95_put_var(nid, nvarid, rain_fall)

    call nf95_redef(nid)
    call nf95_def_var(nid, 'snow_f', nf90_float, idim2, nvarid)
    call nf95_put_att(nid, nvarid, 'title', 'precipitation solide')
    call nf95_enddef(nid)
    call nf95_put_var(nid, nvarid, snow_fall)

    DO nsrf = 1, nbsrf
       IF (nsrf<=99) THEN
          WRITE (str2, '(i2.2)') nsrf
          call nf95_redef(nid)
          call nf95_def_var(nid, 'RUG'//str2, nf90_float, idim2, nvarid)
          call nf95_put_att(nid, nvarid, 'title', &
               'rugosite de surface No.'//str2)
          call nf95_enddef(nid)
       ELSE
          PRINT *, 'Trop de sous-mailles'
          STOP 1
       END IF
       call nf95_put_var(nid, nvarid, frugs(:, nsrf))
    END DO

    DO nsrf = 1, nbsrf
       IF (nsrf<=99) THEN
          WRITE (str2, '(i2.2)') nsrf
          call nf95_redef(nid)
          call nf95_def_var(nid, 'AGESNO'//str2, nf90_float, idim2, nvarid)
          call nf95_put_att(nid, nvarid, 'title', &
               'Age de la neige surface No.'//str2)
          call nf95_enddef(nid)
       ELSE
          PRINT *, 'Trop de sous-mailles'
          STOP 1
       END IF
       call nf95_put_var(nid, nvarid, agesno(:, nsrf))
    END DO

    call nf95_redef(nid)
    call nf95_def_var(nid, 'ZMEA', nf90_float, idim2, nvarid)
    call nf95_enddef(nid)
    call nf95_put_var(nid, nvarid, zmea)

    call nf95_redef(nid)
    call nf95_def_var(nid, 'ZSTD', nf90_float, idim2, nvarid)
    call nf95_enddef(nid)
    call nf95_put_var(nid, nvarid, zstd)
    call nf95_redef(nid)
    call nf95_def_var(nid, 'ZSIG', nf90_float, idim2, nvarid)
    call nf95_enddef(nid)
    call nf95_put_var(nid, nvarid, zsig)
    call nf95_redef(nid)
    call nf95_def_var(nid, 'ZGAM', nf90_float, idim2, nvarid)
    call nf95_enddef(nid)
    call nf95_put_var(nid, nvarid, zgam)
    call nf95_redef(nid)
    call nf95_def_var(nid, 'ZTHE', nf90_float, idim2, nvarid)
    call nf95_enddef(nid)
    call nf95_put_var(nid, nvarid, zthe)
    call nf95_redef(nid)
    call nf95_def_var(nid, 'ZPIC', nf90_float, idim2, nvarid)
    call nf95_enddef(nid)
    call nf95_put_var(nid, nvarid, zpic)
    call nf95_redef(nid)
    call nf95_def_var(nid, 'ZVAL', nf90_float, idim2, nvarid)
    call nf95_enddef(nid)
    call nf95_put_var(nid, nvarid, zval)

    call nf95_redef(nid)
    call nf95_def_var(nid, 'TANCIEN', nf90_float, idim3, nvarid)
    call nf95_enddef(nid)
    call nf95_put_var(nid, nvarid, pack(t_ancien, .true.))

    call nf95_redef(nid)
    call nf95_def_var(nid, 'QANCIEN', nf90_float, idim3, nvarid)
    call nf95_enddef(nid)
    call nf95_put_var(nid, nvarid, pack(q_ancien, .true.))

    call nf95_redef(nid)
    call nf95_def_var(nid, 'RUGMER', nf90_float, idim2, nvarid)
    call nf95_put_att(nid, nvarid, 'title', &
         'Longueur de rugosite sur mer')
    call nf95_enddef(nid)
    call nf95_put_var(nid, nvarid, frugs(:, is_oce))

    call nf95_redef(nid)
    call nf95_def_var(nid, 'CLWCON', nf90_float, idim2, nvarid)
    call nf95_put_att(nid, nvarid, 'title', 'Eau liquide convective')
    call nf95_enddef(nid)
    call nf95_put_var(nid, nvarid, clwcon)

    call nf95_redef(nid)
    call nf95_def_var(nid, 'RNEBCON', nf90_float, idim2, nvarid)
    call nf95_put_att(nid, nvarid, 'title', 'Nebulosite convective')
    call nf95_enddef(nid)
    call nf95_put_var(nid, nvarid, rnebcon)

    call nf95_redef(nid)
    call nf95_def_var(nid, 'RATQS', nf90_float, idim2, nvarid)
    call nf95_put_att(nid, nvarid, 'title', 'Ratqs')
    call nf95_enddef(nid)
    call nf95_put_var(nid, nvarid, ratqs)

    ! run_off_lic_0

    call nf95_redef(nid)
    call nf95_def_var(nid, 'RUNOFFLIC0', nf90_float, idim2, nvarid)
    call nf95_put_att(nid, nvarid, 'title', 'Runofflic0')
    call nf95_enddef(nid)
    call nf95_put_var(nid, nvarid, run_off_lic_0)


    call nf95_close(nid)

  END SUBROUTINE phyredem

end module phyredem_m
