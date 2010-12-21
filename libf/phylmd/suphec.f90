module suphec_m

  implicit none

contains

  SUBROUTINE suphec

    ! From phylmd/suphec.F,v 1.2 2005/06/06 13:16:33

    ! Initialise certaines constantes et parametres physiques.

    use YOMCST, only: rpi, rclum, rhpla, rkbol, rnavo, rday, rea, repsm, &
         rsiyea, rsiday,romega, rg, ra, r1sa, rsigma, r, rmd, rmo3, rmv, rd, &
         rv, rcpd, rcvd, rcpv, rcvv, rkappa, retv, rcw, rcs, rtt, rlvtt, &
         rlstt, rlmlt, ratm, restt, rgamw, rbetw, ralpw, rgams, rbets, ralps, &
         rgamd, rbetd, ralpd
    use yoethf, only: r2es, r3ies, r3les, r4ies, r4les, r5ies, r5les, rhoh2o, &
         rvtmp2

    !------------------------------------------

    PRINT *, 'Call sequence information: suphec'

    ! 1. DEFINE FUNDAMENTAL CONSTANTS

    print *, 'Constants of the ICM'
    RPI=2.*ASIN(1.)
    RCLUM=299792458.
    RHPLA=6.6260755E-34
    RKBOL=1.380658E-23
    RNAVO=6.0221367E+23
    print *, 'Fundamental constants '
    print '(''           PI = '',E13.7,'' -'')', RPI
    print '(''            c = '',E13.7,''m s-1'')', RCLUM
    print '(''            h = '',E13.7,''J s'')', RHPLA
    print '(''            K = '',E13.7,''J K-1'')', RKBOL
    print '(''            N = '',E13.7,''mol-1'')', RNAVO

    ! 2. DEFINE ASTRONOMICAL CONSTANTS

    RDAY=86400.
    REA=149597870000.
    REPSM=0.409093

    RSIYEA=365.25*RDAY*2.*RPI/6.283076
    RSIDAY=RDAY/(1.+RDAY/RSIYEA)
    ROMEGA=2.*RPI/RSIDAY

    print *, 'Astronomical constants '
    print '(''          day = '',E13.7,'' s'')', RDAY
    print '('' half g. axis = '',E13.7,'' m'')', REA
    print '('' mean anomaly = '',E13.7,'' -'')', REPSM
    print '('' sideral year = '',E13.7,'' s'')', RSIYEA
    print '(''  sideral day = '',E13.7,'' s'')', RSIDAY
    print '(''        omega = '',E13.7,'' s-1'')', ROMEGA

    ! 3.    DEFINE GEOIDE.

    RG=9.80665
    RA=6371229.
    R1SA=SNGL(1.D0/DBLE(RA))
    print *, '        Geoide      '
    print '(''      Gravity = '',E13.7,'' m s-2'')', RG
    print '('' Earth radius = '',E13.7,'' m'')', RA
    print '('' Inverse E.R. = '',E13.7,'' m'')', R1SA

    ! 4.    DEFINE RADIATION CONSTANTS.

    rsigma = 2.*rpi**5 * (rkbol/rhpla)**3 * rkbol/rclum/rclum/15.
    print *, '       Radiation    '
    print '('' Stefan-Bol.  = '',E13.7,'' W m-2 K-4'')',   RSIGMA

    ! 5.    DEFINE THERMODYNAMIC CONSTANTS, GAS PHASE.

    R=RNAVO*RKBOL
    RMD=28.9644
    RMO3=47.9942
    RMV=18.0153
    RD=1000.*R/RMD
    RV=1000.*R/RMV
    RCPD=3.5*RD
    RCVD=RCPD-RD
    RCPV=4. *RV
    RCVV=RCPV-RV
    RKAPPA=RD/RCPD
    RETV=RV/RD-1.
    print *, 'Thermodynamic, gas  '
    print '('' Perfect gas  = '',e13.7)',  R
    print '('' Dry air mass = '',e13.7)',  RMD
    print '('' Ozone   mass = '',e13.7)',  RMO3
    print '('' Vapour  mass = '',e13.7)',  RMV
    print '('' Dry air cst. = '',e13.7)',  RD
    print '('' Vapour  cst. = '',e13.7)',  RV
    print '(''         Cpd  = '',e13.7)',  RCPD
    print '(''         Cvd  = '',e13.7)',  RCVD
    print '(''         Cpv  = '',e13.7)',  RCPV
    print '(''         Cvv  = '',e13.7)',  RCVV
    print '(''      Rd/Cpd  = '',e13.7)',  RKAPPA
    print '(''     Rv/Rd-1  = '',e13.7)',  RETV

    ! 6.    DEFINE THERMODYNAMIC CONSTANTS, LIQUID PHASE.

    RCW=RCPV
    print *, 'Thermodynamic, liquid  '
    print '(''         Cw   = '',E13.7)',  RCW

    ! 7.    DEFINE THERMODYNAMIC CONSTANTS, SOLID PHASE.

    RCS=RCPV
    print *, 'thermodynamic, solid'
    print '(''         Cs   = '',E13.7)',  RCS

    ! 8.    DEFINE THERMODYNAMIC CONSTANTS, TRANSITION OF PHASE.

    RTT=273.16
    RLVTT=2.5008E+6
    RLSTT=2.8345E+6
    RLMLT=RLSTT-RLVTT
    RATM=100000.
    print *, 'Thermodynamic, trans.  '
    print '('' Fusion point  = '',E13.7)',  RTT
    print '(''        RLvTt  = '',E13.7)',  RLVTT
    print '(''        RLsTt  = '',E13.7)',  RLSTT
    print '(''        RLMlt  = '',E13.7)',  RLMLT
    print '('' Normal press. = '',E13.7)',  RATM

    ! 9.    SATURATED VAPOUR PRESSURE.

    RESTT=611.14
    RGAMW=(RCW-RCPV)/RV
    RBETW=RLVTT/RV+RGAMW*RTT
    RALPW=LOG(RESTT)+RBETW/RTT+RGAMW*LOG(RTT)
    RGAMS=(RCS-RCPV)/RV
    RBETS=RLSTT/RV+RGAMS*RTT
    RALPS=LOG(RESTT)+RBETS/RTT+RGAMS*LOG(RTT)
    RGAMD=RGAMS-RGAMW
    RBETD=RBETS-RBETW
    RALPD=RALPS-RALPW

    ! calculer les constantes pour les fonctions thermodynamiques

    RVTMP2=RCPV/RCPD-1.
    RHOH2O=RATM/100.
    R2ES=RESTT*RD/RV
    R3LES=17.269
    R3IES=21.875
    R4LES=35.86
    R4IES=7.66
    R5LES=R3LES*(RTT-R4LES)
    R5IES=R3IES*(RTT-R4IES)

  END SUBROUTINE suphec

end module suphec_m
