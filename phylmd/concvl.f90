module concvl_m

  IMPLICIT NONE

contains

  SUBROUTINE concvl(paprs, play, t, q, u, v, sig1, w01, d_t, d_q, d_u, d_v, &
       rain, kbas, itop_con, upwd, dnwd, ma, cape, iflag, qcondc, pmflxr, da, &
       phi, mp)

    ! From phylmd/concvl.F, version 1.3, 2005/04/15 12:36:17
    ! Author: Z. X. Li (LMD/CNRS)
    ! Date: 1993 August 18
    ! Objet : schéma de convection d'Emanuel (1991), interface

    use comconst, only: dtphys
    use cv_driver_m, only: cv_driver
    USE dimphy, ONLY: klev, klon
    USE fcttre, ONLY: foeew
    USE suphec_m, ONLY: retv, rtt
    USE yoethf_m, ONLY: r2es

    REAL, INTENT (IN):: paprs(klon, klev + 1)
    REAL, INTENT (IN):: play(klon, klev)
    REAL, intent(in):: t(klon, klev) ! temperature (K)
    real, intent(in):: q(klon, klev) ! fraction massique de vapeur d'eau
    real, INTENT (IN):: u(klon, klev), v(klon, klev)
    REAL, intent(inout):: sig1(klon, klev), w01(klon, klev)
    REAL, intent(out):: d_t(klon, klev)
    REAL, intent(out):: d_q(klon, klev) ! incr\'ement de la vapeur d'eau
    REAL, intent(out):: d_u(:, :), d_v(:, :) ! (klon, klev)
    REAL, intent(out):: rain(klon) ! pluie (mm / s)
    INTEGER, intent(out):: kbas(klon)
    integer, intent(inout):: itop_con(klon)

    REAL, intent(out):: upwd(klon, klev)
    ! saturated updraft mass flux (kg / m2 / s)

    real, intent(out):: dnwd(klon, klev)
    ! saturated downdraft mass flux (kg / m2 / s)

    REAL ma(klon, klev)
    real cape(klon) ! output (J / kg)
    INTEGER, intent(out):: iflag(klon)
    REAL, intent(out):: qcondc(klon, klev) ! in-cloud water content
    REAL, intent(out):: pmflxr(klon, klev + 1)
    REAL, intent(out):: da(:, :) ! (klon, klev)
    REAL, intent(out):: phi(:, :, :) ! (klon, klev, klev)

    REAL, intent(out):: mp(:, :) ! (klon, klev) Mass flux of the
    ! unsaturated downdraft, defined positive downward, in kg m-2
    ! s-1. M_p in Emanuel (1991 928).

    ! Local:
    REAL zx_qs, cor
    INTEGER i, k
    REAL qs(klon, klev)

    !-----------------------------------------------------------------

    DO k = 1, klev
       DO i = 1, klon
          zx_qs = min(0.5, r2es * foeew(t(i, k), rtt >= t(i, k)) / play(i, k))
          cor = 1. / (1. - retv * zx_qs)
          qs(i, k) = zx_qs * cor
       END DO
    END DO

    CALL cv_driver(t, q, qs, u, v, play / 100., paprs / 100., iflag, d_t, &
         d_q, d_u, d_v, rain, pmflxr, sig1, w01, kbas, itop_con, ma, upwd, &
         dnwd, qcondc, cape, da, phi, mp)
    rain = rain / 86400.
    d_t = dtphys * d_t
    d_q = dtphys * d_q
    d_u = dtphys * d_u
    d_v = dtphys * d_v

  END SUBROUTINE concvl

end module concvl_m
