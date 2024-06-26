module suphec_m

  use jumble, only: pi, twoPI

  use comconst, only: daysec

  implicit none

  ! A1.0 Fundamental constants
  real, parameter, private:: RCLUM = 299792458. ! speed of light, m s-1
  real, parameter, private:: RHPLA = 6.6260755E-34 ! Planck constant, J s
  real, parameter, private:: KBOL = 1.380658E-23 ! Boltzmann constant, in J K-1
  real, parameter, private:: NAVO = 6.0221367E23 ! Avogadro number, in mol-1

  ! A1.1 Astronomical constants

  REAL, parameter, private:: n_sid = 365.25636
  ! Number of days in sideral year. Capderou 2003 k0784, § 4.2.1.

  REAL, parameter:: ROMEGA = twoPI / daysec * (n_sid + 1) / n_sid
  ! Capderou 2003 k0784, equation 4.20

  ! A1.2 Geoide
  real, parameter:: RG = 9.80665 ! acceleration of gravity, in m s-2

  ! A1.3 Radiation
  REAL, parameter:: rsigma = 2. * pi**5 * (kbol / rhpla)**3 * kbol / rclum**2 &
       / 15.

  ! A1.4 Thermodynamic gas phase
  REAL, parameter:: R = NAVO * KBOL ! ideal gas constant, in J K-1 mol-1
  real, parameter:: MV = 18.0153 ! molar mass of water, in g mol-1

  real, parameter:: RV = 1e3 * R / MV
  ! specific ideal gas constant for water vapor, in J K-1 kg-1
  ! (factor 1e3: conversion from g to kg)

  real, parameter:: MD = 28.9644 ! molar mass of dry air, in g mol-1

  real, parameter:: RD = 1e3 * R / MD
  ! specific ideal gas constant for dry air, in J K-1 kg-1
  ! (factor 1e3: conversion from g to kg)

  real, parameter:: RCPV = 4. * RV 
  ! specific heat capacity at constant pressure of water vapor, in J K-1 kg-1

  real, parameter:: RCVV = RCPV - RV
  ! specific heat capacity at constant volume of water vapor, in J K-1 kg-1

  real, parameter:: RCPD = 7. / 2 * RD 
  ! specific heat capacity at constant pressure of dry air, in J K-1 kg-1

  real, parameter:: RCVD = RCPD - RD
  ! specific heat capacity at constant volume of dry air, in J K-1 kg-1
  
  real, parameter:: RMO3 = 47.9942
  REAL, parameter:: RKAPPA = RD / RCPD
  real, parameter:: RETV = RV / RD - 1.

  ! A1.5, 6 Thermodynamic liquid, solid phases

  REAL, parameter:: RCW = RCPV ! LIQUID PHASE Cw

  ! A1.7 Thermodynamic transition of phase
  real, parameter:: RTT = 273.16

  real, parameter:: RLVTT = 2.5008E+6 
  ! specific latent heat of vaporization of water at triple point, in J kg-1

  real, parameter:: RLSTT = 2.8345E+6
  REAL, parameter:: RLMLT = RLSTT - RLVTT
  real, parameter:: RATM = 1e5

  ! A1.8 Curve of saturation
  real, parameter:: RESTT = 611.14

  private pi, twoPI

contains

  SUBROUTINE suphec

    ! From phylmd/suphec.F, version 1.2 2005/06/06 13:16:33

    !------------------------------------------

    PRINT *, 'Call sequence information: suphec'

    print *, 'Astronomical constants '
    print *, 'omega = ', ROMEGA, 's-1'

    print *, 'Radiation constants:'
    print *, 'Stefan-Boltzmann constant = ', RSIGMA, ' W m-2 K-4'

    print *, 'Thermodynamical constants, gas phase:'
    print *, "rd = ", RD, "J K-1 kg-1"
    print *, "rv = ", RV, "J K-1 kg-1"
    print *, 'Cpd = ', RCPD
    print *, 'Cvd = ', RCVD
    print *, 'Cvv = ', RCVV
    print *, 'Rd / Cpd = ', RKAPPA
    print *, 'Rv / Rd - 1 = ', RETV
    print *, 'RCPV = ', RCPV

    print *, 'Thermodynamic, transition of phase:'
    print *, 'RLMlt = ', RLMLT

  END SUBROUTINE suphec

end module suphec_m
