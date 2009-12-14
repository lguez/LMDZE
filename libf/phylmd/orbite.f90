MODULE orbite_m

  ! From phylmd/orbite.F, v 1.1.1.1 2004/05/19 12:53:08

  IMPLICIT NONE

CONTAINS

  SUBROUTINE orbite(xjour, longi, dist)

    USE yomcst, ONLY : r_ecc, r_peri
    use numer_rec, only: pi

    ! Auteur(s): Z.X. Li (LMD/CNRS)
    ! Date: 1993/08/18
    ! Pour un jour donné, calcule la longitude vraie de la Terre (par
    ! rapport au point vernal, 21 mars) dans son orbite solaire. Calcule aussi
    ! la distance Terre-Soleil, c'est-à-dire l'unité astronomique.

    REAL, INTENT (IN):: xjour ! jour de l'année à compter du premier janvier

    REAL, INTENT (OUT):: longi
    ! longitude vraie de la Terre dans son orbite solaire, par
    ! rapport au point vernal (21 mars), en degrés

    REAL, INTENT (OUT), OPTIONAL:: dist
    ! distance terre-soleil (par rapport a la moyenne)

    ! Variables locales
    REAL pir, xl, xllp, xee, xse, xlam, anm, ranm, ranv

    !----------------------------------------------------------------------

    pir = pi / 180.
    xl = r_peri + 180.
    xllp = xl * pir
    xee = r_ecc * r_ecc
    xse = sqrt(1. - xee)
    xlam = (r_ecc / 2 + r_ecc * xee / 8.) * (1. + xse) * sin(xllp) &
         - xee / 4. * (0.5 + xse) * sin(2.*xllp) &
         + r_ecc * xee / 8. * (1. / 3. + xse) * sin(3.*xllp)
    xlam = 2. * xlam / pir
    anm = xlam + (xjour - 81.) - xl
    ranm = anm * pir
    xee = xee * r_ecc
    ranv = ranm + (2. * r_ecc - xee / 4.) * sin(ranm) + &
         5. / 4. * r_ecc * r_ecc * sin(2 * ranm) &
         + 13. / 12. * xee * sin(3.*ranm)

    longi = ranv / pir + xl

    IF (present(dist)) then
       dist = (1 - r_ecc*r_ecc) &
            / (1 + r_ecc*cos(pir*(longi - (r_peri + 180.))))
    end IF

  END SUBROUTINE orbite

  !*****************************************************************

  SUBROUTINE zenang(longi, gmtime, pdtrad, pmu0, frac)

    USE dimphy, ONLY : klon
    USE yomcst, ONLY : r_incl
    USE phyetat0_m, ONLY : rlat, rlon
    USE comconst, ONLY : pi
    use numer_rec, only: assert

    ! Author : O. Boucher (LMD/CNRS), d'après les routines "zenith" et
    ! "angle" de Z.X. Li

    ! Calcule les valeurs moyennes du cos de l'angle zénithal et
    ! l'ensoleillement moyen entre "gmtime1" et "gmtime2" connaissant la
    ! déclinaison, la latitude et la longitude.
    ! Différent de la routine "angle" en ce sens que "zenang" fournit des
    ! moyennes de "pmu0" et non des valeurs instantanées.
    ! Du coup "frac" prend toutes les valeurs entre 0 et 1.

    ! Date : première version le 13 decembre 1994
    !        revu pour  GCM  le 30 septembre 1996

    REAL, INTENT (IN):: longi
    ! (longitude vraie de la terre dans son plan solaire a partir de
    ! l'equinoxe de printemps) (in degrees)

    REAL, INTENT (IN):: gmtime ! temps universel en fraction de jour
    REAL, INTENT (IN):: pdtrad ! pas de temps du rayonnement (secondes)

    REAL, INTENT (OUT):: pmu0(:) ! (klon)
    ! (cosine of mean zenith angle between "gmtime" and "gmtime+pdtrad")

    REAL, INTENT (OUT), OPTIONAL:: frac(:) ! (klon)
    ! (ensoleillement moyen entre gmtime et gmtime+pdtrad)

    ! Variables local to the procedure:

    INTEGER i
    REAL gmtime1, gmtime2
    REAL deux_pi

    REAL omega1, omega2, omega
    ! omega1, omega2 : temps 1 et 2 exprimés en radians avec 0 à midi.
    ! omega : heure en radians du coucher de soleil 
    ! -omega est donc l'heure en radians de lever du soleil

    REAL omegadeb, omegafin
    REAL zfrac1, zfrac2, z1_mu, z2_mu
    REAL lat_sun ! déclinaison en radians
    REAL latr ! latitude du point de grille en radians

    !----------------------------------------------------------------------

    if (present(frac)) call assert((/size(pmu0), size(frac)/) == klon, &
         "zenang")
    
    deux_pi = 2*pi

    lat_sun = asin(sin(longi * pi / 180.) * sin(r_incl * pi / 180.))
    ! Capderou (2003 #784, équation 4.49)

    gmtime1 = gmtime*86400.
    gmtime2 = gmtime*86400. + pdtrad

    DO i = 1, klon
       latr = rlat(i)*pi/180.
       omega = 0.0 !--nuit polaire
       IF (latr>=(pi/2.-lat_sun) .OR. latr<=(-pi/2.-lat_sun)) THEN
          omega = pi ! journee polaire
       END IF
       IF (latr<(pi/2.+lat_sun) .AND. latr>(-pi/2.+lat_sun) .AND. &
            latr<(pi/2.-lat_sun) .AND. latr>(-pi/2.-lat_sun)) THEN
          omega = -tan(latr)*tan(lat_sun)
          omega = acos(omega)
       END IF

       omega1 = gmtime1 + rlon(i)*86400.0/360.0
       omega1 = omega1/86400.0*deux_pi
       omega1 = mod(omega1+deux_pi, deux_pi)
       omega1 = omega1 - pi

       omega2 = gmtime2 + rlon(i)*86400.0/360.0
       omega2 = omega2/86400.0*deux_pi
       omega2 = mod(omega2+deux_pi, deux_pi)
       omega2 = omega2 - pi

       TEST_OMEGA12: IF (omega1<=omega2) THEN
          ! on est dans la meme journee locale
          IF (omega2<=-omega .OR. omega1>=omega .OR. omega<1E-5) THEN
             ! nuit
             IF (present(frac)) frac(i) = 0.0
             pmu0(i) = 0.0
          ELSE
             ! jour + nuit / jour
             omegadeb = max(-omega, omega1)
             omegafin = min(omega, omega2)
             IF (present(frac)) frac(i) = (omegafin-omegadeb)/(omega2-omega1)
             pmu0(i) = sin(latr)*sin(lat_sun) + cos(latr)*cos(lat_sun)*(sin( &
                  omegafin)-sin(omegadeb))/(omegafin-omegadeb)
          END IF
       ELSE TEST_OMEGA12
          !---omega1 GT omega2 -- a cheval sur deux journees
          !-------------------entre omega1 et pi
          IF (omega1>=omega) THEN !--nuit
             zfrac1 = 0.0
             z1_mu = 0.0
          ELSE !--jour+nuit
             omegadeb = max(-omega, omega1)
             omegafin = omega
             zfrac1 = omegafin - omegadeb
             z1_mu = sin(latr)*sin(lat_sun) + cos(latr)*cos(lat_sun)*(sin( &
                  omegafin)-sin(omegadeb))/(omegafin-omegadeb)
          END IF
          !---------------------entre -pi et omega2
          IF (omega2<=-omega) THEN !--nuit
             zfrac2 = 0.0
             z2_mu = 0.0
          ELSE !--jour+nuit
             omegadeb = -omega
             omegafin = min(omega, omega2)
             zfrac2 = omegafin - omegadeb
             z2_mu = sin(latr)*sin(lat_sun) + cos(latr)*cos(lat_sun)*(sin( &
                  omegafin)-sin(omegadeb))/(omegafin-omegadeb)

          END IF
          !-----------------------moyenne 
          IF (present(frac)) frac(i) = (zfrac1+zfrac2)/ &
               (omega2+deux_pi-omega1)
          pmu0(i) = (zfrac1*z1_mu+zfrac2*z2_mu)/max(zfrac1+zfrac2, 1.E-10)
       END IF TEST_OMEGA12
    END DO

  END SUBROUTINE zenang

END MODULE orbite_m
