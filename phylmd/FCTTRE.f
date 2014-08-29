module FCTTRE

  ! From phylmd/FCTTRE.inc, version 1.2, 2004/06/22 11:45:32

  ! This module includes the thermodynamical functions for the cycle
  ! 39 ECMWF physics package. Consistent with "SUPHEC_M" basic
  ! physical constants, assuming the partial pressure of water vapour
  ! is given by a first order Taylor expansion of "Qs(T)" with respect
  ! to temperature, using constants in "yoethf_m".

  ! Probably from Buck, 1981, Journal of Applied Meteorology, volume
  ! 20, number 12, page 1527.

  implicit none

  LOGICAL, PARAMETER:: thermcep = .TRUE.

contains

  elemental REAL function FOEEW(T, DEL)

    use yoethf_m, only: R3LES, R3IES, R4LES, R4IES
    use SUPHEC_M, only: rtt

    REAL, intent(in):: T
    logical, intent(in):: DEL ! ice, else liquid

    !-----------------------

    FOEEW = exp(merge(R3IES / (T - R4IES), R3lES / (T - R4lES), del) &
         * (T - RTT))

  end function FOEEW

  !******************************************

  REAL function FOEDE(T, DEL, P5ARG, QS, PCOARG)

    use yoethf_m, only: R4LES, R4IES

    REAL, intent(in):: T
    logical, intent(in):: DEL ! ice, else liquid
    real, intent(in):: P5ARG, QS, PCOARG

    !-----------------------

    FOEDE = QS * PCOARG * P5ARG / (T - merge(R4IES, R4lES, del))**2

  end function FOEDE

  !******************************************

  elemental REAL function qsats(t)

    REAL, intent(in):: T

    !-----------------------

    qsats = 100. * 0.622 &
         * 10.**(2.07023 - 0.00320991 * t - 2484.896 / t + 3.56654 * LOG10(t))

  end function qsats

  !******************************************

  elemental REAL function qsatl(t)

    REAL, intent(in):: T

    !-----------------------

    qsatl = 100. * 0.622 * 10.**(23.8319 - 2948.964 / t &
         - 5.028 * LOG10(t) - 29810.16 * EXP(- 0.0699382 * t) &
         + 25.21935 * EXP(- 2999.924 / t))

  end function qsatl

  !******************************************

  REAL function dqsats(t, qs)

    use SUPHEC_M, only: RLVTT, rcpd

    REAL, intent(in):: T, qs

    !-----------------------

    dqsats = RLVTT / RCPD * qs * (3.56654/t &
         +2484.896*LOG(10.)/t**2 - 0.00320991*LOG(10.))

  end function dqsats

  !******************************************

  REAL function dqsatl(t, qs)

    use SUPHEC_M, only: RLVTT, rcpd

    REAL, intent(in):: T, qs

    !-----------------------

    dqsatl = RLVTT / RCPD * qs * LOG(10.) &
         * (2948.964/t**2-5.028/LOG(10.)/t &
         +25.21935*2999.924/t**2*EXP(-2999.924/t) &
         +29810.16*0.0699382*EXP(-0.0699382*t))

  end function dqsatl

end module FCTTRE
