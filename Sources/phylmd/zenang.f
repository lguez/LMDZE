module zenang_m

  IMPLICIT NONE

contains

  SUBROUTINE zenang(longi, gmtime, pdtrad, mu0, frac)

    ! Author: O. Boucher (LMD/CNRS), d'après les routines "zenith" et
    ! "angle" de Z.X. Li

    ! Date : première version le 13 décembre 1994, revu pour GCM le 30
    ! septembre 1996

    ! Calcule les valeurs moyennes du cos de l'angle zénithal et
    ! l'ensoleillement moyen entre "gmtime" et "gmtime + pdtrad"
    ! connaissant la déclinaison, la latitude et la longitude.
    ! Différent de la routine "angle" en ce sens que "zenang" fournit
    ! des moyennes de "mu0" et non des valeurs instantanées. Du coup
    ! "frac" prend toutes les valeurs entre 0 et 1. Cf. Capderou (2003
    ! 784, equation 9.11).

    USE dimphy, ONLY: klon
    USE yomcst, ONLY: r_incl
    USE phyetat0_m, ONLY: rlat, rlon
    use nr_util, only: assert, pi, twopi

    REAL, INTENT(IN):: longi
    ! longitude vraie de la terre dans son plan solaire à partir de
    ! l'équinoxe de printemps (in degrees)

    REAL, INTENT(IN):: gmtime ! temps universel en fraction de jour
    REAL, INTENT(IN):: pdtrad ! pas de temps du rayonnement (s)

    REAL, INTENT(OUT):: mu0(:) ! (klon)
    ! cosine of mean zenith angle between "gmtime" and "gmtime+pdtrad"

    REAL, INTENT(OUT), OPTIONAL:: frac(:) ! (klon)
    ! ensoleillement moyen entre gmtime et gmtime+pdtrad

    ! Local:

    INTEGER i
    REAL gmtime1, gmtime2
    REAL omega1, omega2 ! temps 1 et 2 exprimés en radians avec 0 à midi

    REAL omega ! heure en rad du coucher de soleil 
    ! - omega est donc l'heure en rad de lever du soleil

    REAL omegadeb, omegafin
    REAL zfrac1, zfrac2, z1_mu, z2_mu
    REAL lat_sun ! déclinaison en radians
    REAL latr ! latitude du point de grille en radians

    !----------------------------------------------------------------------

    if (present(frac)) call assert((/size(mu0), size(frac)/) == klon, "zenang")

    lat_sun = asin(sin(longi * pi / 180.) * sin(r_incl * pi / 180.))
    ! Capderou (2003 784, equation 4.49)

    gmtime1 = gmtime*86400.
    gmtime2 = gmtime*86400. + pdtrad

    DO i = 1, klon
       latr = rlat(i)*pi/180.
       omega = 0.0 ! nuit polaire
       IF (latr>=(pi/2.-lat_sun) .OR. latr<=(-pi/2.-lat_sun)) THEN
          omega = pi ! journée polaire
       END IF
       IF (latr<(pi/2.+lat_sun) .AND. latr>(-pi/2.+lat_sun) .AND. &
            latr<(pi/2.-lat_sun) .AND. latr>(-pi/2.-lat_sun)) THEN
          omega = -tan(latr)*tan(lat_sun)
          omega = acos(omega)
       END IF

       omega1 = gmtime1 + rlon(i)*86400.0/360.0
       omega1 = omega1/86400.0*twopi
       omega1 = mod(omega1+twopi, twopi)
       omega1 = omega1 - pi

       omega2 = gmtime2 + rlon(i)*86400.0/360.0
       omega2 = omega2/86400.0*twopi
       omega2 = mod(omega2+twopi, twopi)
       omega2 = omega2 - pi

       IF (omega1<=omega2) THEN
          ! on est dans la meme journee locale
          IF (omega2<=-omega .OR. omega1>=omega .OR. omega<1E-5) THEN
             ! nuit
             IF (present(frac)) frac(i) = 0.0
             mu0(i) = 0.0
          ELSE
             ! jour + nuit / jour
             omegadeb = max(-omega, omega1)
             omegafin = min(omega, omega2)
             IF (present(frac)) frac(i) = (omegafin-omegadeb)/(omega2-omega1)
             mu0(i) = sin(latr) * sin(lat_sun) + cos(latr) * cos(lat_sun) &
                  * (sin(omegafin) - sin(omegadeb)) / (omegafin - omegadeb)
          END IF
       ELSE 
          ! omega1 > omega2, à cheval sur deux journées
          ! entre omega1 et pi
          IF (omega1>=omega) THEN ! nuit
             zfrac1 = 0.0
             z1_mu = 0.0
          ELSE ! jour+nuit
             omegadeb = max(-omega, omega1)
             omegafin = omega
             zfrac1 = omegafin - omegadeb
             z1_mu = sin(latr) * sin(lat_sun) + cos(latr) * cos(lat_sun) &
                  * (sin(omegafin) - sin(omegadeb)) / (omegafin - omegadeb)
          END IF
          ! entre -pi et omega2
          IF (omega2<=-omega) THEN ! nuit
             zfrac2 = 0.0
             z2_mu = 0.0
          ELSE ! jour+nuit
             omegadeb = -omega
             omegafin = min(omega, omega2)
             zfrac2 = omegafin - omegadeb
             z2_mu = sin(latr) * sin(lat_sun) + cos(latr) * cos(lat_sun) &
                  * (sin(omegafin) - sin(omegadeb)) / (omegafin - omegadeb)
          END IF
          ! moyenne 
          IF (present(frac)) frac(i) = (zfrac1+zfrac2)/ (omega2+twopi-omega1)
          mu0(i) = (zfrac1*z1_mu+zfrac2*z2_mu)/max(zfrac1+zfrac2, 1.E-10)
       END IF
    END DO

  END SUBROUTINE zenang

end module zenang_m