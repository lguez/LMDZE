module suphec_m

  use nr_util, only: pi, twoPI

  implicit none

  ! A1.0 Fundamental constants
  real, parameter:: RCLUM = 299792458. ! speed of light, m s-1
  real, parameter:: RHPLA = 6.6260755E-34 ! Planck constant, J s
  real, parameter:: KBOL = 1.380658E-23 ! Boltzmann constant, in J K-1
  real, parameter:: NAVO = 6.0221367E23 ! Avogadro number, in mol-1

  ! A1.1 Astronomical constants

  real, parameter:: RDAY = 86400.

  REAL, parameter:: RSIYEA = 365.25 * RDAY * 2. * PI / 6.283076
  ! sideral year, in s

  REAL, parameter:: RSIDAY = RDAY / (1. + RDAY / RSIYEA) ! sideral day, in s
  REAL, parameter:: ROMEGA = twoPI / RSIDAY

  ! A1.2 Geoide
  real, parameter:: RG = 9.80665 ! acceleration of gravity, in m s-2
  real, parameter:: RA = 6371229.

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
  real, save:: RCS

  ! A1.7 Thermodynamic transition of phase
  REAL, save:: RLMLT
  real, parameter:: RTT = 273.16

  real, parameter:: RLVTT = 2.5008E+6 
  ! specific latent heat of vaporization of water at triple point, in J kg-1

  real, parameter:: RLSTT = 2.8345E+6
  real, parameter:: RATM = 1e5

  ! A1.8 Curve of saturation
  REAL, save:: RALPW, RBETW, RGAMW, RALPS, RBETS, RGAMS
  real, parameter:: RESTT = 611.14
  REAL, save:: RALPD, RBETD, RGAMD

  private pi, twoPI

contains

  SUBROUTINE suphec

    ! From phylmd/suphec.F, version 1.2 2005/06/06 13:16:33
    ! Initialise certaines constantes et certains paramètres physiques.

    !------------------------------------------

    PRINT *, 'Call sequence information: suphec'

    print *, 'Astronomical constants '
    print '('' omega = '', E13.7, '' s-1'')', ROMEGA

    print *, 'Geoid:'
    print '('' Gravity = '', E13.7, '' m s-2'')', RG
    print '('' Earth radius = '', E13.7, '' m'')', RA

    print *, 'Radiation constants:'
    print '('' Stefan-Bol. = '', E13.7, '' W m-2 K-4'')', RSIGMA

    print *, 'Thermodynamical constants, gas phase:'
    print *, "rd = ", RD, "J K-1 kg-1"
    print *, "rv = ", RV, "J K-1 kg-1"
    print '('' Cpd = '', e13.7)', RCPD
    print '('' Cvd = '', e13.7)', RCVD
    print '('' Cvv = '', e13.7)', RCVV
    print '('' Rd / Cpd = '', e13.7)', RKAPPA
    print *, 'Rv / Rd - 1 = ', RETV

    ! 7. DEFINE THERMODYNAMIC CONSTANTS, SOLID PHASE.

    RCS = RCPV
    print *, 'thermodynamic, solid'
    print '('' Cs = '', E13.7)', RCS

    ! 8. DEFINE THERMODYNAMIC CONSTANTS, TRANSITION OF PHASE.

    RLMLT = RLSTT-RLVTT
    print *, 'Thermodynamic, transition of phase:'
    print '('' Fusion point = '', E13.7)', RTT
    print '('' RLsTt = '', E13.7)', RLSTT
    print '('' RLMlt = '', E13.7)', RLMLT
    print '('' Normal pressure = '', E13.7)', RATM

    ! 9. SATURATED VAPOUR PRESSURE.

    RGAMW = (RCW-RCPV) / RV
    RBETW = RLVTT / RV+RGAMW*RTT
    RALPW = LOG(RESTT)+RBETW / RTT+RGAMW*LOG(RTT)
    RGAMS = (RCS-RCPV) / RV
    RBETS = RLSTT / RV+RGAMS*RTT
    RALPS = LOG(RESTT)+RBETS / RTT+RGAMS*LOG(RTT)
    RGAMD = RGAMS-RGAMW
    RBETD = RBETS-RBETW
    RALPD = RALPS-RALPW

  END SUBROUTINE suphec

end module suphec_m
