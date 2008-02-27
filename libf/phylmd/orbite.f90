module orbite_m

  ! From phylmd/orbite.F, v 1.1.1.1 2004/05/19 12:53:08

  IMPLICIT none

contains

  SUBROUTINE orbite(xjour, longi, dist)

    use YOMCST

    ! Auteur(s): Z.X. Li (LMD/CNRS) Date: 1993/08/18 Pour un jour
    ! donné, calcule la longitude vraie de la Terre (par rapport au
    ! point vernal, 21 mars) dans son orbite solaire. Calcule aussi
    ! la distance Terre-Soleil, c'est-à-dire l'unité astronomique.

    REAL, intent(in):: xjour ! jour de l'annee a compter du 1er janvier

    real, intent(out):: longi
    ! (longitude vraie de la Terre dans son orbite solaire, par
    ! rapport au point vernal (21 mars), en degrés)

    real, intent(out), optional:: dist
    ! (distance terre-soleil (par rapport a la moyenne))

    ! Variables locales
    REAL pir, xl, xllp, xee, xse, xlam, dlamm, anm, ranm, anv, ranv

    !----------------------------------------------------------------------

    pir = 4.0*ATAN(1.0) / 180.0
    xl=R_peri+180.0
    xllp=xl*pir
    xee=R_ecc*R_ecc
    xse=SQRT(1.0-xee)
    xlam = (R_ecc/2.0+R_ecc*xee/8.0)*(1.0+xse)*SIN(xllp) &
         - xee/4.0*(0.5+xse)*SIN(2.0*xllp) &
         + R_ecc*xee/8.0*(1.0/3.0+xse)*SIN(3.0*xllp)
    xlam=2.0*xlam/pir
    dlamm=xlam+(xjour-81.0)
    anm=dlamm-xl
    ranm=anm*pir
    xee=xee*R_ecc
    ranv=ranm+(2.0*R_ecc-xee/4.0)*SIN(ranm) &
         +5.0/4.0*R_ecc*R_ecc*SIN(2.0*ranm) &
         +13.0/12.0*xee*SIN(3.0*ranm)

    anv=ranv/pir
    longi=anv+xl

    if (present(dist)) dist &
         = (1-R_ecc*R_ecc) /(1+R_ecc*COS(pir*(longi-(R_peri+180.0))))

  END SUBROUTINE orbite

  !*****************************************************************

  SUBROUTINE angle(longi, lati, frac, muzero)

    use dimphy, only: klon
    use YOMCST, only: R_incl

    ! Author: Z.X. Li (LMD/CNRS), date: 1993/08/18
    ! Calcule la durée d'ensoleillement pour un jour et la hauteur du
    ! soleil (cosinus de l'angle zénithal) moyennée sur la journee.

    ! Arguments:
    ! longi----INPUT-R- la longitude vraie de la terre dans son plan 
    !                   solaire a partir de l'equinoxe de printemps (degre)
    ! lati-----INPUT-R- la latitude d'un point sur la terre (degre)
    ! frac-----OUTPUT-R la duree d'ensoleillement dans la journee divisee
    !                   par 24 heures (unite en fraction de 0 a 1)
    ! muzero---OUTPUT-R la moyenne du cosinus de l'angle zinithal sur
    !                   la journee (0 a 1)

    REAL longi
    REAL lati(klon), frac(klon), muzero(klon)
    REAL lat, omega, lon_sun, lat_sun
    REAL pi_local, incl
    INTEGER i

    !----------------------------------------------------------------------

    pi_local = 4.0 * ATAN(1.0)
    incl=R_incl * pi_local / 180.

    lon_sun = longi * pi_local / 180.0
    lat_sun = ASIN (sin(lon_sun)*SIN(incl) )

    DO i = 1, klon
       lat = lati(i) * pi_local / 180.0

       IF ( lat .GE. (pi_local/2.+lat_sun) &
            .OR. lat <= (-pi_local/2.+lat_sun)) THEN
          omega = 0.0   ! nuit polaire
       ELSE IF ( lat.GE.(pi_local/2.-lat_sun) &
            .OR. lat <= (-pi_local/2.-lat_sun)) THEN
          omega = pi_local   ! journee polaire
       ELSE
          omega = -TAN(lat)*TAN(lat_sun)
          omega = ACOS (omega)
       ENDIF

       frac(i) = omega / pi_local

       IF (omega .GT. 0.0) THEN
          muzero(i) = SIN(lat)*SIN(lat_sun) &
               + COS(lat)*COS(lat_sun)*SIN(omega) / omega
       ELSE
          muzero(i) = 0.0
       ENDIF
    ENDDO

  END SUBROUTINE angle

  !*****************************************************************

  SUBROUTINE zenang(longi, gmtime, pdtrad, pmu0, frac)

    use dimphy, only: klon
    use YOMCST, only: r_incl
    use phyetat0_m, only: rlat, rlon

    ! Author : O. Boucher (LMD/CNRS), d'après les routines "zenith" et
    ! "angle" de Z.X. Li

    ! Calcule les valeurs moyennes du cos de l'angle zénithal et
    ! l'ensoleillement moyen entre "gmtime1" et "gmtime2" connaissant la
    ! déclinaison, la latitude et la longitude.
    ! Différent de la routine "angle" en ce sens que "zenang" fournit des
    ! moyennes de "pmu0" et non des valeurs instantanées.
    ! Du coup "frac" prend toutes les valeurs entre 0 et 1.

    ! Date : premiere version le 13 decembre 1994
    !          revu pour  GCM  le 30 septembre 1996

    real, intent(in):: longi
    ! (longitude vraie de la terre dans son plan solaire a partir de
    ! l'equinoxe de printemps (degre))

    real, intent(in):: gmtime ! temps universel en fraction de jour
    real, intent(in):: pdtrad ! pas de temps du rayonnement (secondes)

    real, intent(out):: pmu0(klon)
    ! (cosinus de l'angle zenithal moyen entre gmtime et gmtime+pdtrad)

    real, intent(out), optional:: frac(klon)
    ! (ensoleillement moyen entre gmtime et gmtime+pdtrad)

    ! Variables local to the procedure:

    integer i
    real gmtime1, gmtime2
    real pi_local, deux_pi_local, incl
    real omega1, omega2, omega
    ! omega1, omega2 : temps 1 et 2 exprime en radian avec 0 a midi.
    ! omega : heure en radian du coucher de soleil 
    ! -omega est donc l'heure en radian de lever du soleil
    real omegadeb, omegafin
    real zfrac1, zfrac2, z1_mu, z2_mu
    real lat_sun          ! declinaison en radian
    real lon_sun          ! longitude solaire en radian
    real latr             ! latitude du pt de grille en radian

    !----------------------------------------------------------------------

    pi_local = 4.0 * ATAN(1.0)
    deux_pi_local = 2.0 * pi_local
    incl=R_incl * pi_local / 180.

    lon_sun = longi * pi_local / 180.0
    lat_sun = ASIN (SIN(lon_sun)*SIN(incl) )

    gmtime1=gmtime*86400.
    gmtime2=gmtime*86400.+pdtrad

    DO i = 1, klon
       latr = rlat(i) * pi_local / 180.
       omega=0.0  !--nuit polaire
       IF (latr >= (pi_local/2.-lat_sun) &
            .OR. latr <= (-pi_local/2.-lat_sun)) THEN
          omega = pi_local  ! journee polaire
       ENDIF
       IF (latr < (pi_local/2.+lat_sun).AND. &
            latr.GT.(-pi_local/2.+lat_sun).AND. &
            latr < (pi_local/2.-lat_sun).AND. &
            latr.GT.(-pi_local/2.-lat_sun)) THEN
          omega = -TAN(latr)*TAN(lat_sun)
          omega = ACOS(omega)
       ENDIF

       omega1 = gmtime1 + rlon(i)*86400.0/360.0
       omega1 = omega1 / 86400.0*deux_pi_local
       omega1 = MOD (omega1+deux_pi_local, deux_pi_local)
       omega1 = omega1 - pi_local

       omega2 = gmtime2 + rlon(i)*86400.0/360.0
       omega2 = omega2 / 86400.0*deux_pi_local
       omega2 = MOD (omega2+deux_pi_local, deux_pi_local)
       omega2 = omega2 - pi_local

       test_omega12: IF (omega1 <= omega2) THEN  
          ! on est dans la meme journee locale
          IF (omega2 <= -omega .OR. omega1 >= omega .OR. omega < 1e-5) THEN
             ! nuit
             if (present(frac)) frac(i)=0.0
             pmu0(i)=0.0
          ELSE
             ! jour + nuit / jour
             omegadeb=MAX(-omega, omega1)
             omegafin=MIN(omega, omega2)
             if (present(frac)) frac(i)=(omegafin-omegadeb)/(omega2-omega1)
             pmu0(i)=SIN(latr)*SIN(lat_sun) +  &
                  COS(latr)*COS(lat_sun)* &
                  (SIN(omegafin)-SIN(omegadeb))/ &
                  (omegafin-omegadeb)        
          ENDIF
       ELSE test_omega12
          !---omega1 GT omega2 -- a cheval sur deux journees
          !-------------------entre omega1 et pi
          IF (omega1 >= omega) THEN  !--nuit
             zfrac1=0.0
             z1_mu =0.0
          ELSE                       !--jour+nuit
             omegadeb=MAX(-omega, omega1)
             omegafin=omega
             zfrac1=omegafin-omegadeb
             z1_mu =SIN(latr)*SIN(lat_sun) + &
                  COS(latr)*COS(lat_sun)* &
                  (SIN(omegafin)-SIN(omegadeb))/ &
                  (omegafin-omegadeb)
          ENDIF
          !---------------------entre -pi et omega2
          IF (omega2 <= -omega) THEN   !--nuit
             zfrac2=0.0
             z2_mu =0.0
          ELSE                         !--jour+nuit
             omegadeb=-omega
             omegafin=MIN(omega, omega2)
             zfrac2=omegafin-omegadeb
             z2_mu =SIN(latr)*SIN(lat_sun) + &
                  COS(latr)*COS(lat_sun)* &
                  (SIN(omegafin)-SIN(omegadeb))/ &
                  (omegafin-omegadeb)

          ENDIF
          !-----------------------moyenne 
          if (present(frac)) &
               frac(i)=(zfrac1+zfrac2)/(omega2+deux_pi_local-omega1)
          pmu0(i)=(zfrac1*z1_mu+zfrac2*z2_mu)/MAX(zfrac1+zfrac2, 1.E-10)
       ENDIF test_omega12
    ENDDO

  END SUBROUTINE zenang

  !*****************************************************************

  SUBROUTINE zenith (longi, gmtime, lat, long, pmu0, fract)

    use dimphy, only: klon
    use YOMCST, only: R_incl
    use comconst, only: pi

    ! Author: Z.X. Li (LMD/ENS)

    ! Calcule le cosinus de l'angle zénithal du Soleil en connaissant
    ! la déclinaison du Soleil, la latitude et la longitude du point
    ! sur la Terre, et le temps universel.

    REAL, intent(in):: longi ! déclinaison du soleil (en degrés)

    REAL, intent(in):: gmtime
    ! (temps universel en fraction de jour, entre 0 et 1)

    REAL, intent(in):: lat(klon), long(klon) ! latitude et longitude en degres
    REAL, intent(out):: pmu0(klon) ! cosinus de l'angle zenithal
    REAL, intent(out):: fract(klon)

    ! Variables local to the procedure:

    INTEGER n
    REAL zpir, omega
    REAL lat_sun

    !----------------------------------------------------------------------

    zpir = pi / 180.0
    lat_sun = ASIN (SIN(longi * zpir)*SIN(R_incl * zpir) )

    ! Initialisation à la nuit:
    pmu0 = 0.
    fract = 0.

    ! 1 degree in longitude = 240 s in time

    DO n = 1, klon
       omega = (gmtime + long(n) / 360.) * 2 * pi
       omega = MOD(omega + 2 * pi, 2 * pi)
       omega = omega - pi
       pmu0(n) = sin(lat(n)*zpir) * sin(lat_sun) &
            + cos(lat(n)*zpir) * cos(lat_sun) &
            * cos(omega)
       pmu0(n) = MAX (pmu0(n), 0.)
       IF (pmu0(n) > 1.E-6) fract(n)=1.0
    ENDDO

  END SUBROUTINE zenith

end module orbite_m
