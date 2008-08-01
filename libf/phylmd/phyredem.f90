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
    USE netcdf95, ONLY : nf95_create, nf95_put_att
    USE netcdf, ONLY : nf90_clobber, nf90_global

    INCLUDE 'netcdf.inc'

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
    INTEGER :: ierr

    INTEGER :: isoil, nsrf
    CHARACTER (7) :: str7
    CHARACTER (2) :: str2

    !------------------------------------------------------------

    PRINT *, 'Call sequence information: phyredem'
    CALL nf95_create(fichnom, nf90_clobber, nid)

    call nf95_put_att(nid, nf90_global, 'title', 'Fichier redémarrage physique')
    call nf95_put_att(nid, nf90_global, "itau_phy", itau_phy)

    ierr = nf_def_dim(nid, 'points_physiques', klon, idim2)
    ierr = nf_def_dim(nid, 'horizon_vertical', klon*klev, idim3)

    ierr = nf_def_var(nid, 'longitude', nf_float, 1, idim2, nvarid)
    ierr = nf_put_att_text(nid, nvarid, 'title', 32, &
         'Longitudes de la grille physique')
    ierr = nf_enddef(nid)
    ierr = nf_put_var_real(nid, nvarid, rlon)

    ierr = nf_redef(nid)
    ierr = nf_def_var(nid, 'latitude', nf_float, 1, idim2, nvarid)
    ierr = nf_put_att_text(nid, nvarid, 'title', 31, &
         'Latitudes de la grille physique')
    ierr = nf_enddef(nid)
    ierr = nf_put_var_real(nid, nvarid, rlat)

    ! PB ajout du masque terre/mer

    ierr = nf_redef(nid)
    ierr = nf_def_var(nid, 'masque', nf_float, 1, idim2, nvarid)
    ierr = nf_put_att_text(nid, nvarid, 'title', 16, 'masque terre mer')
    ierr = nf_enddef(nid)
    ierr = nf_put_var_real(nid, nvarid, zmasq)
    ! BP ajout des fraction de chaque sous-surface

    ! 1. fraction de terre

    ierr = nf_redef(nid)
    ierr = nf_def_var(nid, 'FTER', nf_float, 1, idim2, nvarid)
    ierr = nf_put_att_text(nid, nvarid, 'title', 21, 'fraction de continent')
    ierr = nf_enddef(nid)
    ierr = nf_put_var_real(nid, nvarid, pctsrf(1:klon, is_ter))

    ! 2. Fraction de glace de terre

    ierr = nf_redef(nid)
    ierr = nf_def_var(nid, 'FLIC', nf_float, 1, idim2, nvarid)
    ierr = nf_put_att_text(nid, nvarid, 'title', 24, 'fraction glace de terre')
    ierr = nf_enddef(nid)
    ierr = nf_put_var_real(nid, nvarid, pctsrf(1:klon, is_lic))

    ! 3. fraction ocean

    ierr = nf_redef(nid)
    ierr = nf_def_var(nid, 'FOCE', nf_float, 1, idim2, nvarid)
    ierr = nf_put_att_text(nid, nvarid, 'title', 14, 'fraction ocean')
    ierr = nf_enddef(nid)
    ierr = nf_put_var_real(nid, nvarid, pctsrf(1:klon, is_oce))

    ! 4. Fraction glace de mer

    ierr = nf_redef(nid)
    ierr = nf_def_var(nid, 'FSIC', nf_float, 1, idim2, nvarid)
    ierr = nf_put_att_text(nid, nvarid, 'title', 18, 'fraction glace mer')
    ierr = nf_enddef(nid)
    ierr = nf_put_var_real(nid, nvarid, pctsrf(1:klon, is_sic))



    DO nsrf = 1, nbsrf
       IF (nsrf<=99) THEN
          WRITE (str2, '(i2.2)') nsrf
          ierr = nf_redef(nid)
          ierr = nf_def_var(nid, 'TS'//str2, nf_float, 1, idim2, nvarid)
          ierr = nf_put_att_text(nid, nvarid, 'title', 28, &
               'Temperature de surface No.'//str2)
          ierr = nf_enddef(nid)
       ELSE
          PRINT *, 'Trop de sous-mailles'
          STOP 1
       END IF
       ierr = nf_put_var_real(nid, nvarid, tsol(1, nsrf))
    END DO

    DO nsrf = 1, nbsrf
       DO isoil = 1, nsoilmx
          IF (isoil<=99 .AND. nsrf<=99) THEN
             WRITE (str7, '(i2.2, "srf", i2.2)') isoil, nsrf
             ierr = nf_redef(nid)
             ierr = nf_def_var(nid, 'Tsoil'//str7, nf_float, 1, idim2, nvarid)
             ierr = nf_put_att_text(nid, nvarid, 'title', 29, &
                  'Temperature du sol No.'//str7)
             ierr = nf_enddef(nid)
          ELSE
             PRINT *, 'Trop de couches'
             STOP 1
          END IF
          ierr = nf_put_var_real(nid, nvarid, tsoil(1, isoil, nsrf))
       END DO
    END DO

    !IM "slab" ocean
    ierr = nf_redef(nid)
    ierr = nf_def_var(nid, 'TSLAB', nf_float, 1, idim2, nvarid)
    ierr = nf_put_att_text(nid, nvarid, 'title', 33, &
         'Ecart de la SST (pour slab-ocean)')
    ierr = nf_enddef(nid)
    ierr = nf_put_var_real(nid, nvarid, tslab)

    ierr = nf_redef(nid)
    ierr = nf_def_var(nid, 'SEAICE', nf_float, 1, idim2, nvarid)
    ierr = nf_put_att_text(nid, nvarid, 'title', 33, &
         'Glace de mer kg/m2 (pour slab-ocean)')
    ierr = nf_enddef(nid)
    ierr = nf_put_var_real(nid, nvarid, seaice)

    DO nsrf = 1, nbsrf
       IF (nsrf<=99) THEN
          WRITE (str2, '(i2.2)') nsrf
          ierr = nf_redef(nid)
          ierr = nf_def_var(nid, 'QS'//str2, nf_float, 1, idim2, nvarid)
          ierr = nf_put_att_text(nid, nvarid, 'title', 25, &
               'Humidite de surface No.'//str2)
          ierr = nf_enddef(nid)
       ELSE
          PRINT *, 'Trop de sous-mailles'
          STOP 1
       END IF
       ierr = nf_put_var_real(nid, nvarid, qsurf(1, nsrf))
    END DO

    ierr = nf_redef(nid)
    ierr = nf_def_var(nid, 'QSOL', nf_float, 1, idim2, nvarid)
    ierr = nf_put_att_text(nid, nvarid, 'title', 20, 'Eau dans le sol (mm)')
    ierr = nf_enddef(nid)
    ierr = nf_put_var_real(nid, nvarid, qsol)

    DO nsrf = 1, nbsrf
       IF (nsrf<=99) THEN
          WRITE (str2, '(i2.2)') nsrf
          ierr = nf_redef(nid)
          ierr = nf_def_var(nid, 'ALBE'//str2, nf_float, 1, idim2, nvarid)
          ierr = nf_put_att_text(nid, nvarid, 'title', 23, &
               'albedo de surface No.'//str2)
          ierr = nf_enddef(nid)
       ELSE
          PRINT *, 'Trop de sous-mailles'
          STOP 1
       END IF
       ierr = nf_put_var_real(nid, nvarid, albedo(1, nsrf))
    END DO

    !IM BEG albedo LW
    DO nsrf = 1, nbsrf
       IF (nsrf<=99) THEN
          WRITE (str2, '(i2.2)') nsrf
          ierr = nf_redef(nid)
          ierr = nf_def_var(nid, 'ALBLW'//str2, nf_float, 1, idim2, nvarid)
          ierr = nf_put_att_text(nid, nvarid, 'title', 23, &
               'albedo LW de surface No.'//str2)
          ierr = nf_enddef(nid)
       ELSE
          PRINT *, 'Trop de sous-mailles'
          STOP 1
       END IF
       ierr = nf_put_var_real(nid, nvarid, alblw(1, nsrf))
    END DO
    !IM END albedo LW

    DO nsrf = 1, nbsrf
       IF (nsrf<=99) THEN
          WRITE (str2, '(i2.2)') nsrf
          ierr = nf_redef(nid)
          ierr = nf_def_var(nid, 'EVAP'//str2, nf_float, 1, idim2, nvarid)
          ierr = nf_put_att_text(nid, nvarid, 'title', 28, &
               'Evaporation de surface No.'//str2)
          ierr = nf_enddef(nid)
       ELSE
          PRINT *, 'Trop de sous-mailles'
          STOP 1
       END IF
       ierr = nf_put_var_real(nid, nvarid, evap(1, nsrf))
    END DO


    DO nsrf = 1, nbsrf
       IF (nsrf<=99) THEN
          WRITE (str2, '(i2.2)') nsrf
          ierr = nf_redef(nid)
          ierr = nf_def_var(nid, 'SNOW'//str2, nf_float, 1, idim2, nvarid)
          ierr = nf_put_att_text(nid, nvarid, 'title', 22, &
               'Neige de surface No.'//str2)
          ierr = nf_enddef(nid)
       ELSE
          PRINT *, 'Trop de sous-mailles'
          STOP 1
       END IF
       ierr = nf_put_var_real(nid, nvarid, snow(1, nsrf))
    END DO


    ierr = nf_redef(nid)
    ierr = nf_def_var(nid, 'RADS', nf_float, 1, idim2, nvarid)
    ierr = nf_put_att_text(nid, nvarid, 'title', 28, &
         'Rayonnement net a la surface')
    ierr = nf_enddef(nid)
    ierr = nf_put_var_real(nid, nvarid, radsol)

    ierr = nf_redef(nid)
    ierr = nf_def_var(nid, 'solsw', nf_float, 1, idim2, nvarid)
    ierr = nf_put_att_text(nid, nvarid, 'title', 32, &
         'Rayonnement solaire a la surface')
    ierr = nf_enddef(nid)
    ierr = nf_put_var_real(nid, nvarid, solsw)

    ierr = nf_redef(nid)
    ierr = nf_def_var(nid, 'sollw', nf_float, 1, idim2, nvarid)
    ierr = nf_put_att_text(nid, nvarid, 'title', 27, &
         'Rayonnement IF a la surface')
    ierr = nf_enddef(nid)
    ierr = nf_put_var_real(nid, nvarid, sollw)

    ierr = nf_redef(nid)
    ierr = nf_def_var(nid, 'fder', nf_float, 1, idim2, nvarid)
    ierr = nf_put_att_text(nid, nvarid, 'title', 14, 'Derive de flux')
    ierr = nf_enddef(nid)
    ierr = nf_put_var_real(nid, nvarid, fder)

    ierr = nf_redef(nid)
    ierr = nf_def_var(nid, 'rain_f', nf_float, 1, idim2, nvarid)
    ierr = nf_put_att_text(nid, nvarid, 'title', 21, 'precipitation liquide')
    ierr = nf_enddef(nid)
    ierr = nf_put_var_real(nid, nvarid, rain_fall)

    ierr = nf_redef(nid)
    ierr = nf_def_var(nid, 'snow_f', nf_float, 1, idim2, nvarid)
    ierr = nf_put_att_text(nid, nvarid, 'title', 20, 'precipitation solide')
    ierr = nf_enddef(nid)
    ierr = nf_put_var_real(nid, nvarid, snow_fall)

    DO nsrf = 1, nbsrf
       IF (nsrf<=99) THEN
          WRITE (str2, '(i2.2)') nsrf
          ierr = nf_redef(nid)
          ierr = nf_def_var(nid, 'RUG'//str2, nf_float, 1, idim2, nvarid)
          ierr = nf_put_att_text(nid, nvarid, 'title', 23, &
               'rugosite de surface No.'//str2)
          ierr = nf_enddef(nid)
       ELSE
          PRINT *, 'Trop de sous-mailles'
          STOP 1
       END IF
       ierr = nf_put_var_real(nid, nvarid, frugs(1, nsrf))
    END DO

    DO nsrf = 1, nbsrf
       IF (nsrf<=99) THEN
          WRITE (str2, '(i2.2)') nsrf
          ierr = nf_redef(nid)
          ierr = nf_def_var(nid, 'AGESNO'//str2, nf_float, 1, idim2, nvarid)
          ierr = nf_put_att_text(nid, nvarid, 'title', 15, &
               'Age de la neige surface No.'//str2)
          ierr = nf_enddef(nid)
       ELSE
          PRINT *, 'Trop de sous-mailles'
          STOP 1
       END IF
       ierr = nf_put_var_real(nid, nvarid, agesno(1, nsrf))
    END DO

    ierr = nf_redef(nid)
    ierr = nf_def_var(nid, 'ZMEA', nf_float, 1, idim2, nvarid)
    ierr = nf_enddef(nid)
    ierr = nf_put_var_real(nid, nvarid, zmea)

    ierr = nf_redef(nid)
    ierr = nf_def_var(nid, 'ZSTD', nf_float, 1, idim2, nvarid)
    ierr = nf_enddef(nid)
    ierr = nf_put_var_real(nid, nvarid, zstd)
    ierr = nf_redef(nid)
    ierr = nf_def_var(nid, 'ZSIG', nf_float, 1, idim2, nvarid)
    ierr = nf_enddef(nid)
    ierr = nf_put_var_real(nid, nvarid, zsig)
    ierr = nf_redef(nid)
    ierr = nf_def_var(nid, 'ZGAM', nf_float, 1, idim2, nvarid)
    ierr = nf_enddef(nid)
    ierr = nf_put_var_real(nid, nvarid, zgam)
    ierr = nf_redef(nid)
    ierr = nf_def_var(nid, 'ZTHE', nf_float, 1, idim2, nvarid)
    ierr = nf_enddef(nid)
    ierr = nf_put_var_real(nid, nvarid, zthe)
    ierr = nf_redef(nid)
    ierr = nf_def_var(nid, 'ZPIC', nf_float, 1, idim2, nvarid)
    ierr = nf_enddef(nid)
    ierr = nf_put_var_real(nid, nvarid, zpic)
    ierr = nf_redef(nid)
    ierr = nf_def_var(nid, 'ZVAL', nf_float, 1, idim2, nvarid)
    ierr = nf_enddef(nid)
    ierr = nf_put_var_real(nid, nvarid, zval)

    ierr = nf_redef(nid)
    ierr = nf_def_var(nid, 'TANCIEN', nf_float, 1, idim3, nvarid)
    ierr = nf_enddef(nid)
    ierr = nf_put_var_real(nid, nvarid, t_ancien)

    ierr = nf_redef(nid)
    ierr = nf_def_var(nid, 'QANCIEN', nf_float, 1, idim3, nvarid)
    ierr = nf_enddef(nid)
    ierr = nf_put_var_real(nid, nvarid, q_ancien)

    ierr = nf_redef(nid)
    ierr = nf_def_var(nid, 'RUGMER', nf_float, 1, idim2, nvarid)
    ierr = nf_put_att_text(nid, nvarid, 'title', 28, &
         'Longueur de rugosite sur mer')
    ierr = nf_enddef(nid)
    ierr = nf_put_var_real(nid, nvarid, frugs(1, is_oce))

    ierr = nf_redef(nid)
    ierr = nf_def_var(nid, 'CLWCON', nf_float, 1, idim2, nvarid)
    ierr = nf_put_att_text(nid, nvarid, 'title', 28, 'Eau liquide convective')
    ierr = nf_enddef(nid)
    ierr = nf_put_var_real(nid, nvarid, clwcon)

    ierr = nf_redef(nid)
    ierr = nf_def_var(nid, 'RNEBCON', nf_float, 1, idim2, nvarid)
    ierr = nf_put_att_text(nid, nvarid, 'title', 28, 'Nebulosite convective')
    ierr = nf_enddef(nid)
    ierr = nf_put_var_real(nid, nvarid, rnebcon)

    ierr = nf_redef(nid)
    ierr = nf_def_var(nid, 'RATQS', nf_float, 1, idim2, nvarid)
    ierr = nf_put_att_text(nid, nvarid, 'title', 5, 'Ratqs')
    ierr = nf_enddef(nid)
    ierr = nf_put_var_real(nid, nvarid, ratqs)

    ! run_off_lic_0

    ierr = nf_redef(nid)
    ierr = nf_def_var(nid, 'RUNOFFLIC0', nf_float, 1, idim2, nvarid)
    ierr = nf_put_att_text(nid, nvarid, 'title', 10, 'Runofflic0')
    ierr = nf_enddef(nid)
    ierr = nf_put_var_real(nid, nvarid, run_off_lic_0)


    ierr = nf_close(nid)

  END SUBROUTINE phyredem

end module phyredem_m
