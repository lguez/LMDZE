cIM ctes ds clesphys.h   SUBROUTINE SW(PSCT, RCO2, PRMU0, PFRAC, 
      SUBROUTINE SW(PSCT, PRMU0, PFRAC, 
     S              PPMB, PDP, 
     S              PPSOL, PALBD, PALBP,
     S              PTAVE, PWV, PQS, POZON, PAER,
     S              PCLDSW, PTAU, POMEGA, PCG,
     S              PHEAT, PHEAT0,
     S              PALBPLA,PTOPSW,PSOLSW,PTOPSW0,PSOLSW0,
     S              ZFSUP,ZFSDN,ZFSUP0,ZFSDN0,
     S              tauae, pizae, cgae,
     s              PTAUA, POMEGAA,
     S              PTOPSWAD,PSOLSWAD,PTOPSWAI,PSOLSWAI,
     J              ok_ade, ok_aie )
      
      use dimens_m
      use dimphy
      use clesphys
      use YOMCST
      use raddim
      IMPLICIT none

C
C     ------------------------------------------------------------------
C
C     PURPOSE.
C     --------
C
C          THIS ROUTINE COMPUTES THE SHORTWAVE RADIATION FLUXES IN TWO
C     SPECTRAL INTERVALS FOLLOWING FOUQUART AND BONNEL (1980).
C
C     METHOD.
C     -------
C
C          1. COMPUTES ABSORBER AMOUNTS                 (SWU)
C          2. COMPUTES FLUXES IN 1ST SPECTRAL INTERVAL  (SW1S)
C          3. COMPUTES FLUXES IN 2ND SPECTRAL INTERVAL  (SW2S)
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
C        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 89-07-14
C        95-01-01   J.-J. MORCRETTE  Direct/Diffuse Albedo
c        03-11-27   J. QUAAS Introduce aerosol forcings (based on BOUCHER)
C     ------------------------------------------------------------------
C
C* ARGUMENTS:
C
      REAL*8 PSCT  ! constante solaire (valeur conseillee: 1370)
cIM ctes ds clesphys.h   REAL*8 RCO2  ! concentration CO2 (IPCC: 353.E-06*44.011/28.97)
C
      REAL*8 PPSOL(KDLON)        ! SURFACE PRESSURE (PA)
      REAL*8 PDP(KDLON,KFLEV)    ! LAYER THICKNESS (PA)
      REAL*8 PPMB(KDLON,KFLEV+1) ! HALF-LEVEL PRESSURE (MB)
C
      REAL*8 PRMU0(KDLON)  ! COSINE OF ZENITHAL ANGLE
      REAL*8 PFRAC(KDLON)  ! fraction de la journee
C
      REAL*8 PTAVE(KDLON,KFLEV)  ! LAYER TEMPERATURE (K)
      REAL*8 PWV(KDLON,KFLEV)    ! SPECIFIC HUMIDITY (KG/KG)
      REAL*8 PQS(KDLON,KFLEV)    ! SATURATED WATER VAPOUR (KG/KG)
      REAL*8 POZON(KDLON,KFLEV)  ! OZONE CONCENTRATION (KG/KG)
      REAL*8 PAER(KDLON,KFLEV,5) ! AEROSOLS' OPTICAL THICKNESS
C
      REAL*8 PALBD(KDLON,2)  ! albedo du sol (lumiere diffuse)
      REAL*8 PALBP(KDLON,2)  ! albedo du sol (lumiere parallele)
C
      REAL*8 PCLDSW(KDLON,KFLEV)    ! CLOUD FRACTION
      REAL*8 PTAU(KDLON,2,KFLEV)    ! CLOUD OPTICAL THICKNESS
      REAL*8 PCG(KDLON,2,KFLEV)     ! ASYMETRY FACTOR
      REAL*8 POMEGA(KDLON,2,KFLEV)  ! SINGLE SCATTERING ALBEDO
C
      REAL*8 PHEAT(KDLON,KFLEV) ! SHORTWAVE HEATING (K/DAY)
      REAL*8 PHEAT0(KDLON,KFLEV)! SHORTWAVE HEATING (K/DAY) clear-sky
      REAL*8 PALBPLA(KDLON)     ! PLANETARY ALBEDO
      REAL*8 PTOPSW(KDLON)      ! SHORTWAVE FLUX AT T.O.A.
      REAL*8 PSOLSW(KDLON)      ! SHORTWAVE FLUX AT SURFACE
      REAL*8 PTOPSW0(KDLON)     ! SHORTWAVE FLUX AT T.O.A. (CLEAR-SKY)
      REAL*8 PSOLSW0(KDLON)     ! SHORTWAVE FLUX AT SURFACE (CLEAR-SKY)
C
C* LOCAL VARIABLES:
C
      REAL*8 ZOZ(KDLON,KFLEV)
      REAL*8 ZAKI(KDLON,2)     
      REAL*8 ZCLD(KDLON,KFLEV)
      REAL*8 ZCLEAR(KDLON) 
      REAL*8 ZDSIG(KDLON,KFLEV)
      REAL*8 ZFACT(KDLON)
      REAL*8 ZFD(KDLON,KFLEV+1)
      REAL*8 ZFDOWN(KDLON,KFLEV+1)
      REAL*8 ZFU(KDLON,KFLEV+1)
      REAL*8 ZFUP(KDLON,KFLEV+1)
      REAL*8 ZRMU(KDLON)
      REAL*8 ZSEC(KDLON)
      REAL*8 ZUD(KDLON,5,KFLEV+1)
      REAL*8 ZCLDSW0(KDLON,KFLEV)
c
      REAL*8 ZFSUP(KDLON,KFLEV+1)
      REAL*8 ZFSDN(KDLON,KFLEV+1)
      REAL*8 ZFSUP0(KDLON,KFLEV+1)
      REAL*8 ZFSDN0(KDLON,KFLEV+1)
C
      INTEGER inu, jl, jk, i, k, kpl1
c
      INTEGER swpas  ! Every swpas steps, sw is calculated
      PARAMETER(swpas=1)
c
      INTEGER itapsw
      LOGICAL appel1er
      DATA itapsw /0/
      DATA appel1er /.TRUE./
cjq-Introduced for aerosol forcings
      real*8 flag_aer
      logical ok_ade, ok_aie    ! use aerosol forcings or not?
      real*8 tauae(kdlon,kflev,2)  ! aerosol optical properties
      real*8 pizae(kdlon,kflev,2)  ! (see aeropt.F)
      real*8 cgae(kdlon,kflev,2)   ! -"-
      REAL*8 PTAUA(KDLON,2,KFLEV)    ! CLOUD OPTICAL THICKNESS (pre-industrial value)
      REAL*8 POMEGAA(KDLON,2,KFLEV)  ! SINGLE SCATTERING ALBEDO
      REAL*8 PTOPSWAD(KDLON)     ! SHORTWAVE FLUX AT T.O.A.(+AEROSOL DIR)
      REAL*8 PSOLSWAD(KDLON)     ! SHORTWAVE FLUX AT SURFACE(+AEROSOL DIR)
      REAL*8 PTOPSWAI(KDLON)     ! SHORTWAVE FLUX AT T.O.A.(+AEROSOL IND)
      REAL*8 PSOLSWAI(KDLON)     ! SHORTWAVE FLUX AT SURFACE(+AEROSOL IND)
cjq - Fluxes including aerosol effects
      REAL*8 ZFSUPAD(KDLON,KFLEV+1)
      REAL*8 ZFSDNAD(KDLON,KFLEV+1)
      REAL*8 ZFSUPAI(KDLON,KFLEV+1)
      REAL*8 ZFSDNAI(KDLON,KFLEV+1)
      logical initialized
      SAVE ZFSUPAD, ZFSDNAD, ZFSUPAI, ZFSDNAI ! aerosol fluxes
!rv
      save flag_aer
      data initialized/.false./
cjq-end
      if(.not.initialized) then
        flag_aer=0.
        initialized=.TRUE.
      endif
!rv
      
c
      IF (appel1er) THEN
         PRINT*, 'SW calling frequency : ', swpas
         PRINT*, "   In general, it should be 1"
         appel1er = .FALSE.
      ENDIF
C     ------------------------------------------------------------------
      IF (MOD(itapsw,swpas).EQ.0) THEN
c
      DO JK = 1 , KFLEV
      DO JL = 1, KDLON
         ZCLDSW0(JL,JK) = 0.0
         IF (bug_ozone) then
           ZOZ(JL,JK) = POZON(JL,JK)*46.6968/RG
     .               *PDP(JL,JK)*(101325.0/PPSOL(JL))
         ELSE
c        Correction MPL 100505
           ZOZ(JL,JK) = POZON(JL,JK)*RMD/RMO3*46.6968/RG*PDP(JL,JK)
         ENDIF           
      ENDDO
      ENDDO
C
C
c clear-sky:
cIM ctes ds clesphys.h  CALL SWU(PSCT,RCO2,ZCLDSW0,PPMB,PPSOL,
      CALL SWU(PSCT,ZCLDSW0,PPMB,PPSOL,
     S         PRMU0,PFRAC,PTAVE,PWV,
     S         ZAKI,ZCLD,ZCLEAR,ZDSIG,ZFACT,ZRMU,ZSEC,ZUD)
      INU = 1
      CALL SW1S(INU,
     S     PAER, flag_aer, tauae, pizae, cgae,
     S     PALBD, PALBP, PCG, ZCLD, ZCLEAR, ZCLDSW0,
     S     ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,
     S     ZFD, ZFU)
      INU = 2
      CALL SW2S(INU,
     S     PAER, flag_aer, tauae, pizae, cgae,
     S     ZAKI, PALBD, PALBP, PCG, ZCLD, ZCLEAR, ZCLDSW0,
     S     ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,
     S     PWV, PQS,
     S     ZFDOWN, ZFUP)
      DO JK = 1 , KFLEV+1
      DO JL = 1, KDLON
         ZFSUP0(JL,JK) = (ZFUP(JL,JK)   + ZFU(JL,JK)) * ZFACT(JL)
         ZFSDN0(JL,JK) = (ZFDOWN(JL,JK) + ZFD(JL,JK)) * ZFACT(JL)
      ENDDO
      ENDDO
      
      flag_aer=0.0
      CALL SWU(PSCT,PCLDSW,PPMB,PPSOL,
     S         PRMU0,PFRAC,PTAVE,PWV,
     S         ZAKI,ZCLD,ZCLEAR,ZDSIG,ZFACT,ZRMU,ZSEC,ZUD)
      INU = 1
      CALL SW1S(INU,
     S     PAER, flag_aer, tauae, pizae, cgae,
     S     PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,
     S     ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,
     S     ZFD, ZFU)
      INU = 2
      CALL SW2S(INU,
     S     PAER, flag_aer, tauae, pizae, cgae,
     S     ZAKI, PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,
     S     ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,
     S     PWV, PQS,
     S    ZFDOWN, ZFUP)

c cloudy-sky:
      
      DO JK = 1 , KFLEV+1
      DO JL = 1, KDLON
         ZFSUP(JL,JK) = (ZFUP(JL,JK)   + ZFU(JL,JK)) * ZFACT(JL)
         ZFSDN(JL,JK) = (ZFDOWN(JL,JK) + ZFD(JL,JK)) * ZFACT(JL)
      ENDDO
      ENDDO
      
c      
      IF (ok_ade) THEN
c
c cloudy-sky + aerosol dir OB
      flag_aer=1.0
      CALL SWU(PSCT,PCLDSW,PPMB,PPSOL,
     S         PRMU0,PFRAC,PTAVE,PWV,
     S         ZAKI,ZCLD,ZCLEAR,ZDSIG,ZFACT,ZRMU,ZSEC,ZUD)
      INU = 1
      CALL SW1S(INU,
     S     PAER, flag_aer, tauae, pizae, cgae,
     S     PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,
     S     ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,
     S     ZFD, ZFU)
      INU = 2
      CALL SW2S(INU,
     S     PAER, flag_aer, tauae, pizae, cgae,
     S     ZAKI, PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,
     S     ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,
     S     PWV, PQS,
     S    ZFDOWN, ZFUP)
      DO JK = 1 , KFLEV+1
      DO JL = 1, KDLON
         ZFSUPAD(JL,JK) = ZFSUP(JL,JK) 
         ZFSDNAD(JL,JK) = ZFSDN(JL,JK) 
         ZFSUP(JL,JK) = (ZFUP(JL,JK)   + ZFU(JL,JK)) * ZFACT(JL)
         ZFSDN(JL,JK) = (ZFDOWN(JL,JK) + ZFD(JL,JK)) * ZFACT(JL)
      ENDDO
      ENDDO 
      
      ENDIF ! ok_ade
      
      IF (ok_aie) THEN
         
cjq   cloudy-sky + aerosol direct + aerosol indirect
      flag_aer=1.0
      CALL SWU(PSCT,PCLDSW,PPMB,PPSOL,
     S         PRMU0,PFRAC,PTAVE,PWV,
     S         ZAKI,ZCLD,ZCLEAR,ZDSIG,ZFACT,ZRMU,ZSEC,ZUD)
      INU = 1
      CALL SW1S(INU,
     S     PAER, flag_aer, tauae, pizae, cgae,
     S     PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,
     S     ZDSIG, POMEGAA, ZOZ, ZRMU, ZSEC, PTAUA, ZUD,
     S     ZFD, ZFU)
      INU = 2
      CALL SW2S(INU,
     S     PAER, flag_aer, tauae, pizae, cgae,
     S     ZAKI, PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,
     S     ZDSIG, POMEGAA, ZOZ, ZRMU, ZSEC, PTAUA, ZUD,
     S     PWV, PQS,
     S    ZFDOWN, ZFUP)
      DO JK = 1 , KFLEV+1
      DO JL = 1, KDLON
         ZFSUPAI(JL,JK) = ZFSUP(JL,JK) 
         ZFSDNAI(JL,JK) = ZFSDN(JL,JK)          
         ZFSUP(JL,JK) = (ZFUP(JL,JK)   + ZFU(JL,JK)) * ZFACT(JL)
         ZFSDN(JL,JK) = (ZFDOWN(JL,JK) + ZFD(JL,JK)) * ZFACT(JL)
      ENDDO
      ENDDO
      ENDIF ! ok_aie      
cjq -end
      
      itapsw = 0
      ENDIF
      itapsw = itapsw + 1
C
      DO k = 1, KFLEV
         kpl1 = k+1
         DO i = 1, KDLON
            PHEAT(i,k) = -(ZFSUP(i,kpl1)-ZFSUP(i,k))
     .                     -(ZFSDN(i,k)-ZFSDN(i,kpl1))
            PHEAT(i,k) = PHEAT(i,k) * RDAY*RG/RCPD / PDP(i,k)
            PHEAT0(i,k) = -(ZFSUP0(i,kpl1)-ZFSUP0(i,k))
     .                     -(ZFSDN0(i,k)-ZFSDN0(i,kpl1))
            PHEAT0(i,k) = PHEAT0(i,k) * RDAY*RG/RCPD / PDP(i,k)
         ENDDO
      ENDDO
      DO i = 1, KDLON
         PALBPLA(i) = ZFSUP(i,KFLEV+1)/(ZFSDN(i,KFLEV+1)+1.0e-20)
c
         PSOLSW(i) = ZFSDN(i,1) - ZFSUP(i,1)
         PTOPSW(i) = ZFSDN(i,KFLEV+1) - ZFSUP(i,KFLEV+1)
c
         PSOLSW0(i) = ZFSDN0(i,1) - ZFSUP0(i,1)
         PTOPSW0(i) = ZFSDN0(i,KFLEV+1) - ZFSUP0(i,KFLEV+1)
c-OB
         PSOLSWAD(i) = ZFSDNAD(i,1) - ZFSUPAD(i,1)
         PTOPSWAD(i) = ZFSDNAD(i,KFLEV+1) - ZFSUPAD(i,KFLEV+1)
c
         PSOLSWAI(i) = ZFSDNAI(i,1) - ZFSUPAI(i,1)
         PTOPSWAI(i) = ZFSDNAI(i,KFLEV+1) - ZFSUPAI(i,KFLEV+1)
c-fin 
      ENDDO
C
      RETURN
      END
