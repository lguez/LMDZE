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

contains

  elemental REAL function FOEEW(T, ICE)

    use yoethf_m, only: R3LES, R3IES, R4LES, R4IES
    use SUPHEC_M, only: rtt

    REAL, intent(in):: T
    logical, intent(in):: ICE ! else liquid

    !-----------------------

    FOEEW = exp(merge(R3IES / (T - R4IES), R3lES / (T - R4lES), ice) &
         * (T - RTT))

  end function FOEEW

  !******************************************

  REAL function FOEDE(T, ICE, P5ARG, QS, PCOARG)

    use yoethf_m, only: R4LES, R4IES

    REAL, intent(in):: T
    logical, intent(in):: ICE ! else liquid
    real, intent(in):: P5ARG, QS, PCOARG

    !-----------------------

    FOEDE = QS * PCOARG * P5ARG / (T - merge(R4IES, R4lES, ice))**2

  end function FOEDE

end module FCTTRE
