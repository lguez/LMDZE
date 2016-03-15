module concvl_m

  IMPLICIT NONE

contains

  SUBROUTINE concvl(dtime, paprs, play, t, q, u, v, sig1, w01, d_t, d_q, d_u, &
       d_v, rain, snow_con, kbas, itop_con, upwd, dnwd, dnwd0, ma, cape, &
       iflag, qcondc, wd, pmflxr, da, phi, mp)

    ! From phylmd/concvl.F, version 1.3, 2005/04/15 12:36:17
    ! Author: Z. X. Li (LMD/CNRS)
    ! Date: 1993 August 18
    ! Objet : schéma de convection d'Emanuel (1991), interface
    ! (driver commun aux versions 3 et 4)

    use cv_driver_m, only: cv_driver
    USE dimphy, ONLY: klev, klon
    USE fcttre, ONLY: foeew
    USE suphec_m, ONLY: retv, rtt
    USE yoethf_m, ONLY: r2es

    REAL, INTENT (IN):: dtime ! pas d'integration (s)
    REAL, INTENT (IN):: paprs(klon, klev + 1)
    REAL, INTENT (IN):: play(klon, klev)
    REAL, intent(in):: t(klon, klev)
    real, intent(in):: q(klon, klev) ! vapeur d'eau (en kg / kg)
    real, INTENT (IN):: u(klon, klev), v(klon, klev)
    REAL, intent(inout):: sig1(klon, klev), w01(klon, klev)
    REAL, intent(out):: d_t(klon, klev)
    REAL, intent(out):: d_q(klon, klev) ! increment de la vapeur d'eau
    REAL, intent(out):: d_u(klon, klev), d_v(klon, klev)
    REAL, intent(out):: rain(klon) ! pluie (mm / s)
    REAL, intent(out):: snow_con(klon) ! neige (mm / s)
    INTEGER, intent(out):: kbas(klon)
    integer itop_con(klon)

    REAL, intent(out):: upwd(klon, klev)
    ! saturated updraft mass flux (kg / m2 / s)

    real, intent(out):: dnwd(klon, klev)
    ! saturated downdraft mass flux (kg / m2 / s)

    real, intent(out):: dnwd0(klon, klev)
    ! unsaturated downdraft mass flux (kg / m2 / s)

    REAL ma(klon, klev)
    real cape(klon) ! output (J / kg)
    INTEGER iflag(klon)
    REAL qcondc(klon, klev)
    REAL wd(klon)
    REAL pmflxr(klon, klev + 1)
    REAL, intent(inout):: da(klon, klev), phi(klon, klev, klev), mp(klon, klev)

    ! Local:
    REAL zx_qs, cor
    INTEGER i, k
    REAL qs(klon, klev)
    REAL, save:: cbmf(klon)
    INTEGER:: ifrst = 0

    !-----------------------------------------------------------------

    snow_con = 0.

    IF (ifrst == 0) THEN
       ifrst = 1
       cbmf = 0.
    END IF

    DO k = 1, klev
       DO i = 1, klon
          zx_qs = min(0.5, r2es * foeew(t(i, k), rtt >= t(i, k)) / play(i, k))
          cor = 1. / (1. - retv * zx_qs)
          qs(i, k) = zx_qs * cor
       END DO
    END DO

    CALL cv_driver(t, q, qs, u, v, play / 100., paprs / 100., iflag, d_t, &
         d_q, d_u, d_v, rain, pmflxr, cbmf, sig1, w01, kbas, itop_con, dtime, &
         ma, upwd, dnwd, dnwd0, qcondc, wd, cape, da, phi, mp)

    rain = rain / 86400.
    d_t = dtime * d_t
    d_q = dtime * d_q
    d_u = dtime * d_u
    d_v = dtime * d_v

  END SUBROUTINE concvl

end module concvl_m
