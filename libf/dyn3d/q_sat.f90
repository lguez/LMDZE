module q_sat_m

  IMPLICIT none

contains

  elemental real function q_sat(temp, pres)

    ! From dyn3d/q_sat.F, version 1.1.1.1 2004/05/19 12:53:05

    ! Authors : Z. X. Li (LMD/CNRS).
    ! Réécriture vectorisée par F. Hourdin.
    ! This procedure computes the mass fraction of saturating water
    ! vapor, with the formula of ECMWF.

    REAL, intent(in):: temp ! temperature, in K
    REAL, intent(in):: pres ! total pressure, in hPa

    ! Variables local to the procedure:

    REAL, PARAMETER:: r2es = 611.14 * 18.0153 / 28.9644

    REAL r3es
    REAL, PARAMETER:: R3LES=17.269
    REAL, PARAMETER:: R3IES=21.875

    REAL r4es
    REAL, PARAMETER:: R4LES=35.86
    REAL, PARAMETER:: R4IES=7.66

    REAL, PARAMETER:: rtt = 273.16
    REAL, PARAMETER:: retv = 28.9644 / 18.0153 - 1.

    !------------------------------------------------------------------

    IF (temp <= rtt) THEN
       r3es = r3ies
       r4es = r4ies
    ELSE
       r3es = r3les
       r4es = r4les
    ENDIF

    q_sat = r2es / pres * EXP(r3es * (temp - rtt) / (temp - r4es))
    q_sat = MIN(0.5, q_sat)
    q_sat = q_sat / (1. - retv * q_sat)

  END function q_sat

end module q_sat_m
