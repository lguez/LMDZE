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
    use yoethf

    LOGICAL:: firstcall = .TRUE.

    !------------------------------------------

    IF (firstcall) THEN
       PRINT *, 'suphec initialise les constantes du GCM'
       firstcall = .FALSE.

       !*       1.    DEFINE FUNDAMENTAL CONSTANTS.

       WRITE(UNIT=6,FMT='(''0*** Constants of the ICM   ***'')')
       RPI=2.*ASIN(1.)
       RCLUM=299792458.
       RHPLA=6.6260755E-34
       RKBOL=1.380658E-23
       RNAVO=6.0221367E+23
       WRITE(UNIT=6,FMT='('' *** Fundamental constants ***'')')
       WRITE(UNIT=6,FMT='(''           PI = '',E13.7,'' -'')')RPI
       WRITE(UNIT=6,FMT='(''            c = '',E13.7,''m s-1'')') &
            RCLUM
       WRITE(UNIT=6,FMT='(''            h = '',E13.7,''J s'')') &
            RHPLA
       WRITE(UNIT=6,FMT='(''            K = '',E13.7,''J K-1'')') &
            RKBOL
       WRITE(UNIT=6,FMT='(''            N = '',E13.7,''mol-1'')') &
            RNAVO

       !*       2.    DEFINE ASTRONOMICAL CONSTANTS.

       RDAY=86400.
       REA=149597870000.
       REPSM=0.409093

       RSIYEA=365.25*RDAY*2.*RPI/6.283076
       RSIDAY=RDAY/(1.+RDAY/RSIYEA)
       ROMEGA=2.*RPI/RSIDAY

       WRITE(UNIT=6,FMT='('' *** Astronomical constants ***'')')
       WRITE(UNIT=6,FMT='(''          day = '',E13.7,'' s'')')RDAY
       WRITE(UNIT=6,FMT='('' half g. axis = '',E13.7,'' m'')')REA
       WRITE(UNIT=6,FMT='('' mean anomaly = '',E13.7,'' -'')')REPSM
       WRITE(UNIT=6,FMT='('' sideral year = '',E13.7,'' s'')')RSIYEA
       WRITE(UNIT=6,FMT='(''  sideral day = '',E13.7,'' s'')')RSIDAY
       WRITE(UNIT=6,FMT='(''        omega = '',E13.7,'' s-1'')') &
            ROMEGA

       !*       3.    DEFINE GEOIDE.

       RG=9.80665
       RA=6371229.
       R1SA=SNGL(1.D0/DBLE(RA))
       WRITE(UNIT=6,FMT='('' ***         Geoide         ***'')')
       WRITE(UNIT=6,FMT='(''      Gravity = '',E13.7,'' m s-2'')') &
            RG
       WRITE(UNIT=6,FMT='('' Earth radius = '',E13.7,'' m'')')RA
       WRITE(UNIT=6,FMT='('' Inverse E.R. = '',E13.7,'' m'')')R1SA

       !*       4.    DEFINE RADIATION CONSTANTS.

       rsigma = 2.*rpi**5 * (rkbol/rhpla)**3 * rkbol/rclum/rclum/15.
       !IM init. dans conf_phys.F90   RI0=1365.
       WRITE(UNIT=6,FMT='('' ***        Radiation       ***'')')
       WRITE(UNIT=6,FMT='('' Stefan-Bol.  = '',E13.7,'' W m-2 K-4'')')  RSIGMA

       !*       5.    DEFINE THERMODYNAMIC CONSTANTS, GAS PHASE.

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
       WRITE(UNIT=6,FMT='('' *** Thermodynamic, gas     ***'')')
       WRITE(UNIT=6,FMT='('' Perfect gas  = '',e13.7)') R
       WRITE(UNIT=6,FMT='('' Dry air mass = '',e13.7)') RMD
       WRITE(UNIT=6,FMT='('' Ozone   mass = '',e13.7)') RMO3
       WRITE(UNIT=6,FMT='('' Vapour  mass = '',e13.7)') RMV
       WRITE(UNIT=6,FMT='('' Dry air cst. = '',e13.7)') RD
       WRITE(UNIT=6,FMT='('' Vapour  cst. = '',e13.7)') RV
       WRITE(UNIT=6,FMT='(''         Cpd  = '',e13.7)') RCPD
       WRITE(UNIT=6,FMT='(''         Cvd  = '',e13.7)') RCVD
       WRITE(UNIT=6,FMT='(''         Cpv  = '',e13.7)') RCPV
       WRITE(UNIT=6,FMT='(''         Cvv  = '',e13.7)') RCVV
       WRITE(UNIT=6,FMT='(''      Rd/Cpd  = '',e13.7)') RKAPPA
       WRITE(UNIT=6,FMT='(''     Rv/Rd-1  = '',e13.7)') RETV

       !*       6.    DEFINE THERMODYNAMIC CONSTANTS, LIQUID PHASE.

       RCW=RCPV
       WRITE(UNIT=6,FMT='('' *** Thermodynamic, liquid  ***'')')
       WRITE(UNIT=6,FMT='(''         Cw   = '',E13.7)') RCW

       !*       7.    DEFINE THERMODYNAMIC CONSTANTS, SOLID PHASE.

       RCS=RCPV
       WRITE(UNIT=6,FMT='('' *** thermodynamic, solid   ***'')')
       WRITE(UNIT=6,FMT='(''         Cs   = '',E13.7)') RCS

       !*       8.    DEFINE THERMODYNAMIC CONSTANTS, TRANSITION OF PHASE.

       RTT=273.16
       RLVTT=2.5008E+6
       RLSTT=2.8345E+6
       RLMLT=RLSTT-RLVTT
       RATM=100000.
       WRITE(UNIT=6,FMT='('' *** Thermodynamic, trans.  ***'')')
       WRITE(UNIT=6,FMT='('' Fusion point  = '',E13.7)') RTT
       WRITE(UNIT=6,FMT='(''        RLvTt  = '',E13.7)') RLVTT
       WRITE(UNIT=6,FMT='(''        RLsTt  = '',E13.7)') RLSTT
       WRITE(UNIT=6,FMT='(''        RLMlt  = '',E13.7)') RLMLT
       WRITE(UNIT=6,FMT='('' Normal press. = '',E13.7)') RATM
       WRITE(UNIT=6,FMT='('' Latent heat :  '')')

       !*       9.    SATURATED VAPOUR PRESSURE.

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
    ELSE
       PRINT *, 'suphec DEJA APPELE '
    ENDIF

  END SUBROUTINE suphec

end module suphec_m
