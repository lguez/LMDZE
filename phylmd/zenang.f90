module zenang_m

  IMPLICIT NONE

contains

  SUBROUTINE zenang(longi, gmtime, pdtrad, mu0, fract)

    ! Author: O. Boucher (LMD/CNRS), d'après les routines "zenith" et
    ! "angle" de Z. X. Li

    ! Date : première version le 13 décembre 1994, revu pour GCM le 30
    ! septembre 1996

    ! Calcule les valeurs moyennes du cosinus de l'angle zénithal et
    ! l'ensoleillement moyen entre "gmtime" et "gmtime + pdtrad"
    ! connaissant la déclinaison, la latitude et la longitude.
    ! Différent de la routine "angle" parce que "zenang" fournit des
    ! moyennes de "mu0" et non des valeurs instantanées. Du coup
    ! "fract" prend toutes les valeurs entre 0 et 1. Cf. Capderou (2003
    ! k0784, equation 9.11).

    use jumble, only: assert, pi, twopi

    USE dimphy, ONLY: klon
    USE yomcst, ONLY: incl
    USE phyetat0_m, ONLY: rlat, rlon

    REAL, INTENT(IN):: longi
    ! longitude vraie de la terre dans son plan solaire à partir de
    ! l'équinoxe de printemps (in degrees)

    REAL, INTENT(IN):: gmtime ! temps universel en fraction de jour
    REAL, INTENT(IN):: pdtrad ! pas de temps du rayonnement (s)

    REAL, INTENT(OUT):: mu0(:) ! (klon)
    ! cosine of mean zenith angle between "gmtime" and "gmtime + pdtrad"

    REAL, INTENT(OUT), OPTIONAL:: fract(:) ! (klon)
    ! ensoleillement moyen entre gmtime et gmtime + pdtrad

    ! Local:

    INTEGER i
    REAL gmtime1, gmtime2
    REAL omega1, omega2 ! temps 1 et 2 exprimés en radians avec 0 à midi

    REAL omega ! heure en rad du coucher de soleil 
    ! "- omega" est donc l'heure en rad de lever du soleil.

    REAL omegadeb, omegafin
    REAL zfrac1, zfrac2, z1_mu, z2_mu
    REAL lat_sun ! déclinaison en radians
    REAL latr ! latitude du point de grille en radians

    !----------------------------------------------------------------------

    if (present(fract)) call assert((/size(mu0), size(fract)/) == klon, &
         "zenang")

    lat_sun = asin(sin(longi * pi / 180.) * sin(incl * pi / 180.))
    ! Capderou (2003 k0784, equation 4.49)

    gmtime1 = gmtime * 86400.
    gmtime2 = gmtime1 + pdtrad

    DO i = 1, klon
       latr = rlat(i) * pi / 180.
       
       IF (latr >= pi / 2. - lat_sun .OR. latr <= - pi / 2. - lat_sun) then
          omega = pi ! journée polaire
       else IF (latr < pi / 2. + lat_sun .AND. latr > - pi / 2. + lat_sun) THEN
          omega = acos(- tan(latr) * tan(lat_sun))
       else
          omega = 0. ! nuit polaire
       END IF

       omega1 = mod((gmtime1 + rlon(i) * 86400. / 360.) / 86400. * twopi &
            + twopi, twopi) - pi
       omega2 = mod((gmtime2 + rlon(i) * 86400. / 360.) / 86400. * twopi &
            + twopi, twopi) - pi

       IF (omega1 <= omega2) THEN
          ! on est dans la meme journee locale
          IF (omega2 <= - omega .OR. omega1 >= omega .OR. omega < 1E-5) THEN
             ! nuit
             IF (present(fract)) fract(i) = 0.
             mu0(i) = 0.
          ELSE
             ! jour + nuit / jour
             omegadeb = max(- omega, omega1)
             omegafin = min(omega, omega2)
             IF (present(fract)) fract(i) = (omegafin - omegadeb) &
                  / (omega2 - omega1)
             mu0(i) = sin(latr) * sin(lat_sun) + cos(latr) * cos(lat_sun) &
                  * (sin(omegafin) - sin(omegadeb)) / (omegafin - omegadeb)
          END IF
       ELSE 
          ! omega1 > omega2, à cheval sur deux journées
          ! entre omega1 et pi
          IF (omega1 >= omega) THEN ! nuit
             zfrac1 = 0.
             z1_mu = 0.
          ELSE ! jour + nuit
             omegadeb = max(- omega, omega1)
             omegafin = omega
             zfrac1 = omegafin - omegadeb
             z1_mu = sin(latr) * sin(lat_sun) + cos(latr) * cos(lat_sun) &
                  * (sin(omegafin) - sin(omegadeb)) / (omegafin - omegadeb)
          END IF
          ! entre - pi et omega2
          IF (omega2 <= - omega) THEN ! nuit
             zfrac2 = 0.
             z2_mu = 0.
          ELSE ! jour + nuit
             omegadeb = - omega
             omegafin = min(omega, omega2)
             zfrac2 = omegafin - omegadeb
             z2_mu = sin(latr) * sin(lat_sun) + cos(latr) * cos(lat_sun) &
                  * (sin(omegafin) - sin(omegadeb)) / (omegafin - omegadeb)
          END IF
          ! moyenne 
          IF (present(fract)) fract(i) = (zfrac1 + zfrac2) &
               / (omega2 + twopi - omega1)
          mu0(i) = (zfrac1 * z1_mu + zfrac2 * z2_mu) &
               / max(zfrac1 + zfrac2, 1e-10)
       END IF
    END DO

  END SUBROUTINE zenang

end module zenang_m
