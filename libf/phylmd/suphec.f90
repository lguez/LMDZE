module suphec_m

  implicit none

  ! A1.0 Fundamental constants
  REAL RPI
  real, parameter:: RCLUM = 299792458.
  real, parameter:: RHPLA = 6.6260755E-34
  real, parameter:: RKBOL = 1.380658E-23 ! Boltzmann constant, in J K-1
  real, parameter:: RNAVO = 6.0221367E+23 ! Avogadro number, in mol-1

  ! A1.1 Astronomical constants
  REAL RSIYEA, RSIDAY, ROMEGA
  real, parameter:: RDAY = 86400.
  real, parameter:: REA = 149597870000.
  real, parameter:: REPSM = 0.409093

  ! A1.2 Geoide
  real, parameter:: RG = 9.80665 ! acceleration of gravity, in m s-2
  real, parameter:: RA = 6371229.

  ! A1.3 Radiation
  REAL RSIGMA

  ! A1.4 Thermodynamic gas phase
  REAL, parameter:: R = RNAVO * RKBOL ! ideal gas constant, in J K-1 mol-1
  real RV, RCPD, RCPV, RCVD, RCVV
  real, parameter:: RMD = 28.9644 ! molar mass of dry air, in g mol-1

  real, parameter:: RD = 1000. * R / RMD
  ! specific ideal gas constant for dry air, in J K-1 kg-1
  ! (factor 1000: conversion from g to kg)

  real, parameter:: RMO3 = 47.9942
  real, parameter:: RMV = 18.0153
  REAL RKAPPA, RETV

  ! A1.5, 6 Thermodynamic liquid, solid phases
  REAL RCW, RCS

  ! A1.7 Thermodynamic transition of phase
  REAL RLMLT
  real, parameter:: RTT = 273.16
  real, parameter:: RLVTT = 2.5008E+6
  real, parameter:: RLSTT = 2.8345E+6
  real, parameter:: RATM = 100000.

  ! A1.8 Curve of saturation
  REAL RALPW, RBETW, RGAMW, RALPS, RBETS, RGAMS
  real, parameter:: RESTT = 611.14
  REAL RALPD, RBETD, RGAMD

  save

contains

  SUBROUTINE suphec

    ! From phylmd/suphec.F, v 1.2 2005/06/06 13:16:33
    ! Initialise certaines constantes et parametres physiques.

    !------------------------------------------

    PRINT *, 'Call sequence information: suphec'

    ! 1. DEFINE FUNDAMENTAL CONSTANTS

    print *, 'Constants of the ICM'
    RPI = 2.*ASIN(1.)
    print *, 'Fundamental constants '
    print '('' PI = '', E13.7, '' -'')', RPI
    print '('' c = '', E13.7, ''m s-1'')', RCLUM
    print '('' h = '', E13.7, ''J s'')', RHPLA
    print '('' K = '', E13.7, ''J K-1'')', RKBOL
    print '('' N = '', E13.7, ''mol-1'')', RNAVO

    ! 2. DEFINE ASTRONOMICAL CONSTANTS

    RSIYEA = 365.25*RDAY*2.*RPI/6.283076
    RSIDAY = RDAY/(1.+RDAY/RSIYEA)
    ROMEGA = 2.*RPI/RSIDAY

    print *, 'Astronomical constants '
    print '('' day = '', E13.7, '' s'')', RDAY
    print '('' half g. axis = '', E13.7, '' m'')', REA
    print '('' mean anomaly = '', E13.7, '' -'')', REPSM
    print '('' sideral year = '', E13.7, '' s'')', RSIYEA
    print '('' sideral day = '', E13.7, '' s'')', RSIDAY
    print '('' omega = '', E13.7, '' s-1'')', ROMEGA

    ! 3. DEFINE GEOIDE.

    print *, ' Geoide '
    print '('' Gravity = '', E13.7, '' m s-2'')', RG
    print '('' Earth radius = '', E13.7, '' m'')', RA

    ! 4. DEFINE RADIATION CONSTANTS.

    rsigma = 2.*rpi**5 * (rkbol/rhpla)**3 * rkbol/rclum/rclum/15.
    print *, ' Radiation '
    print '('' Stefan-Bol. = '', E13.7, '' W m-2 K-4'')', RSIGMA

    ! 5. DEFINE THERMODYNAMIC CONSTANTS, GAS PHASE.

    RV = 1000.*R/RMV
    RCPD = 3.5*RD
    RCVD = RCPD-RD
    RCPV = 4. *RV
    RCVV = RCPV-RV
    RKAPPA = RD/RCPD
    RETV = RV/RD-1.
    print *, 'Thermodynamic, gas '
    print '('' Perfect gas = '', e13.7)', R
    print '('' Ozone mass = '', e13.7)', RMO3
    print '('' Vapour mass = '', e13.7)', RMV
    print '('' Dry air constant = '', e13.7)', RD
    print '('' Vapour constant = '', e13.7)', RV
    print '('' Cpd = '', e13.7)', RCPD
    print '('' Cvd = '', e13.7)', RCVD
    print '('' Cpv = '', e13.7)', RCPV
    print '('' Cvv = '', e13.7)', RCVV
    print '('' Rd/Cpd = '', e13.7)', RKAPPA
    print '('' Rv/Rd-1 = '', e13.7)', RETV

    ! 6. DEFINE THERMODYNAMIC CONSTANTS, LIQUID PHASE.

    RCW = RCPV
    print *, 'Thermodynamic, liquid '
    print '('' Cw = '', E13.7)', RCW

    ! 7. DEFINE THERMODYNAMIC CONSTANTS, SOLID PHASE.

    RCS = RCPV
    print *, 'thermodynamic, solid'
    print '('' Cs = '', E13.7)', RCS

    ! 8. DEFINE THERMODYNAMIC CONSTANTS, TRANSITION OF PHASE.

    RLMLT = RLSTT-RLVTT
    print *, 'Thermodynamic, trans. '
    print '('' Fusion point = '', E13.7)', RTT
    print '('' RLvTt = '', E13.7)', RLVTT
    print '('' RLsTt = '', E13.7)', RLSTT
    print '('' RLMlt = '', E13.7)', RLMLT
    print '('' Normal press. = '', E13.7)', RATM

    ! 9. SATURATED VAPOUR PRESSURE.

    RGAMW = (RCW-RCPV)/RV
    RBETW = RLVTT/RV+RGAMW*RTT
    RALPW = LOG(RESTT)+RBETW/RTT+RGAMW*LOG(RTT)
    RGAMS = (RCS-RCPV)/RV
    RBETS = RLSTT/RV+RGAMS*RTT
    RALPS = LOG(RESTT)+RBETS/RTT+RGAMS*LOG(RTT)
    RGAMD = RGAMS-RGAMW
    RBETD = RBETS-RBETW
    RALPD = RALPS-RALPW

  END SUBROUTINE suphec

end module suphec_m
