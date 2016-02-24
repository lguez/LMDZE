module concvl_m

  IMPLICIT NONE

contains

  SUBROUTINE concvl(dtime, paprs, play, t, q, u, v, sig1, w01, d_t, d_q, d_u, &
       d_v, rain, snow, kbas, ktop, upwd, dnwd, dnwd0, ma, cape, iflag, &
       qcondc, wd, pmflxr, pmflxs, da, phi, mp)

    ! From phylmd/concvl.F, version 1.3 2005/04/15 12:36:17
    ! Author: Z. X. Li (LMD/CNRS)
    ! Date: 1993 August 18
    ! Objet : schéma de convection d'Emanuel (1991), interface
    ! (driver commun aux versions 3 et 4)

    use clesphys2, only: iflag_con
    use cv_driver_m, only: cv_driver
    USE dimens_m, ONLY: nqmx
    USE dimphy, ONLY: klev, klon
    USE fcttre, ONLY: foeew
    USE suphec_m, ONLY: retv, rtt
    USE yoethf_m, ONLY: r2es

    REAL, INTENT (IN):: dtime ! pas d'integration (s)
    REAL, INTENT (IN):: paprs(klon, klev+1)
    REAL, INTENT (IN):: play(klon, klev)
    REAL, intent(in):: t(klon, klev)
    real, intent(in):: q(klon, klev) ! vapeur d'eau (en kg/kg)
    real, INTENT (IN):: u(klon, klev), v(klon, klev)
    REAL, intent(inout):: sig1(klon, klev), w01(klon, klev)
    REAL, intent(out):: d_t(klon, klev)
    REAL, intent(out):: d_q(klon, klev) ! increment de la vapeur d'eau
    REAL, intent(out):: d_u(klon, klev), d_v(klon, klev)
    REAL, intent(out):: rain(klon) ! pluie (mm/s)
    REAL, intent(out):: snow(klon) ! neige (mm/s)
    INTEGER kbas(klon), ktop(klon)

    REAL, intent(out):: upwd(klon, klev)
    ! saturated updraft mass flux (kg/m**2/s)

    real, intent(out):: dnwd(klon, klev)
    ! saturated downdraft mass flux (kg/m**2/s)

    real, intent(out):: dnwd0(klon, klev)
    ! unsaturated downdraft mass flux (kg/m**2/s)

    REAL ma(klon, klev), cape(klon)
    ! Cape----output-R-CAPE (J/kg)

    INTEGER iflag(klon)
    REAL qcondc(klon, klev)
    REAL wd(klon)
    REAL pmflxr(klon, klev+1), pmflxs(klon, klev+1)
    REAL, intent(inout):: da(klon, klev), phi(klon, klev, klev), mp(klon, klev)

    ! Local:

    REAL em_ph(klon, klev+1), em_p(klon, klev)
    REAL zx_t, zx_qs, zcor
    INTEGER i, k
    REAL qs(klon, klev)
    REAL, save:: cbmf(klon)
    INTEGER:: ifrst = 0

    !-----------------------------------------------------------------

    snow = 0

    IF (ifrst==0) THEN
       ifrst = 1
       DO i = 1, klon
          cbmf(i) = 0.
       END DO
    END IF

    DO k = 1, klev + 1
       DO i = 1, klon
          em_ph(i, k) = paprs(i, k)/100.0
          pmflxs(i, k) = 0.
       END DO
    END DO

    DO k = 1, klev
       DO i = 1, klon
          em_p(i, k) = play(i, k)/100.0
       END DO
    END DO


    IF (iflag_con==4) THEN
       DO k = 1, klev
          DO i = 1, klon
             zx_t = t(i, k)
             zx_qs = min(0.5, r2es*foeew(zx_t, rtt >= zx_t)/em_p(i, k)/100.0)
             zcor = 1./(1.-retv*zx_qs)
             qs(i, k) = zx_qs*zcor
          END DO
       END DO
    ELSE
       ! iflag_con=3 (modification de puristes qui fait la
       ! diff\'erence pour la convergence numerique)
       DO k = 1, klev
          DO i = 1, klon
             zx_t = t(i, k)
             zx_qs = r2es*foeew(zx_t, rtt >= zx_t)/em_p(i, k)/100.0
             zx_qs = min(0.5, zx_qs)
             zcor = 1./(1.-retv*zx_qs)
             zx_qs = zx_qs*zcor
             qs(i, k) = zx_qs
          END DO
       END DO
    END IF

    CALL cv_driver(t, q, qs, u, v, em_p, em_ph, iflag, d_t, d_q, d_u, d_v, &
         rain, pmflxr, cbmf, sig1, w01, kbas, ktop, dtime, ma, upwd, dnwd, &
         dnwd0, qcondc, wd, cape, da, phi, mp)

    DO i = 1, klon
       rain(i) = rain(i)/86400.
    END DO

    DO k = 1, klev
       DO i = 1, klon
          d_t(i, k) = dtime*d_t(i, k)
          d_q(i, k) = dtime*d_q(i, k)
          d_u(i, k) = dtime*d_u(i, k)
          d_v(i, k) = dtime*d_v(i, k)
       END DO
    END DO

  END SUBROUTINE concvl

end module concvl_m
