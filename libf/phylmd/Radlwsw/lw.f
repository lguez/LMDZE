cIM ctes ds clesphys.h   SUBROUTINE LW(RCO2,RCH4,RN2O,RCFC11,RCFC12,
      SUBROUTINE LW(
     .              PPMB, PDP,
     .              PPSOL,PDT0,PEMIS,
     .              PTL, PTAVE, PWV, POZON, PAER,
     .              PCLDLD,PCLDLU,
     .              PVIEW,
     .              PCOLR, PCOLR0,
     .              PTOPLW,PSOLLW,PTOPLW0,PSOLLW0,
     .              psollwdown,
     .              plwup, plwdn, plwup0, plwdn0)
      use dimens_m
      use dimphy
      use clesphys
      use SUPHEC_M
      use raddim
            use raddimlw
      IMPLICIT none
C
C-----------------------------------------------------------------------
C     METHOD.
C     -------
C
C          1. COMPUTES THE PRESSURE AND TEMPERATURE WEIGHTED AMOUNTS OF
C     ABSORBERS.
C          2. COMPUTES THE PLANCK FUNCTIONS ON THE INTERFACES AND THE
C     GRADIENT OF PLANCK FUNCTIONS IN THE LAYERS.
C          3. PERFORMS THE VERTICAL INTEGRATION DISTINGUISHING THE CON-
C     TRIBUTIONS OF THE ADJACENT AND DISTANT LAYERS AND THOSE FROM THE
C     BOUNDARIES.
C          4. COMPUTES THE CLEAR-SKY DOWNWARD AND UPWARD EMISSIVITIES.
C          5. INTRODUCES THE EFFECTS OF THE CLOUDS ON THE FLUXES.
C
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 89-07-14
C-----------------------------------------------------------------------
cIM ctes ds clesphys.h
c     REAL*8 RCO2   ! CO2 CONCENTRATION (IPCC:353.E-06* 44.011/28.97)
c     REAL*8 RCH4   ! CH4 CONCENTRATION (IPCC: 1.72E-06* 16.043/28.97)
c     REAL*8 RN2O   ! N2O CONCENTRATION (IPCC: 310.E-09* 44.013/28.97)
c     REAL*8 RCFC11 ! CFC11 CONCENTRATION (IPCC: 280.E-12* 137.3686/28.97)
c     REAL*8 RCFC12 ! CFC12 CONCENTRATION (IPCC: 484.E-12* 120.9140/28.97)
      REAL*8 PCLDLD(KDLON,KFLEV)  ! DOWNWARD EFFECTIVE CLOUD COVER
      REAL*8 PCLDLU(KDLON,KFLEV)  ! UPWARD EFFECTIVE CLOUD COVER
      REAL*8 PDP(KDLON,KFLEV)     ! LAYER PRESSURE THICKNESS (Pa)
      REAL*8 PDT0(KDLON)          ! SURFACE TEMPERATURE DISCONTINUITY (K)
      REAL*8 PEMIS(KDLON)         ! SURFACE EMISSIVITY
      REAL*8 PPMB(KDLON,KFLEV+1)  ! HALF LEVEL PRESSURE (mb)
      REAL*8 PPSOL(KDLON)         ! SURFACE PRESSURE (Pa)
      REAL*8 POZON(KDLON,KFLEV)   ! O3 CONCENTRATION (kg/kg)
      REAL*8 PTL(KDLON,KFLEV+1)   ! HALF LEVEL TEMPERATURE (K)
      REAL*8 PAER(KDLON,KFLEV,5)  ! OPTICAL THICKNESS OF THE AEROSOLS
      REAL*8 PTAVE(KDLON,KFLEV)   ! LAYER TEMPERATURE (K)
      REAL*8 PVIEW(KDLON)         ! COSECANT OF VIEWING ANGLE
      REAL*8 PWV(KDLON,KFLEV)     ! SPECIFIC HUMIDITY (kg/kg)
C
      REAL*8 PCOLR(KDLON,KFLEV)   ! LONG-WAVE TENDENCY (K/day)
      REAL*8 PCOLR0(KDLON,KFLEV)  ! LONG-WAVE TENDENCY (K/day) clear-sky
      REAL*8 PTOPLW(KDLON)        ! LONGWAVE FLUX AT T.O.A.
      REAL*8 PSOLLW(KDLON)        ! LONGWAVE FLUX AT SURFACE
      REAL*8 PTOPLW0(KDLON)       ! LONGWAVE FLUX AT T.O.A. (CLEAR-SKY)
      REAL*8 PSOLLW0(KDLON)       ! LONGWAVE FLUX AT SURFACE (CLEAR-SKY)
c Rajout LF
      real*8 psollwdown(kdlon)    ! LONGWAVE downwards flux at surface
cIM
      REAL*8 plwup(KDLON,KFLEV+1)  ! LW up total sky
      REAL*8 plwup0(KDLON,KFLEV+1) ! LW up clear sky
      REAL*8 plwdn(KDLON,KFLEV+1)  ! LW down total sky
      REAL*8 plwdn0(KDLON,KFLEV+1) ! LW down clear sky
C-------------------------------------------------------------------------
      REAL*8 ZABCU(KDLON,NUA,3*KFLEV+1)
      REAL*8 ZOZ(KDLON,KFLEV)
c
      REAL*8 ZFLUX(KDLON,2,KFLEV+1) ! RADIATIVE FLUXES (1:up; 2:down)
      REAL*8 ZFLUC(KDLON,2,KFLEV+1) ! CLEAR-SKY RADIATIVE FLUXES
      REAL*8 ZBINT(KDLON,KFLEV+1)            ! Intermediate variable
      REAL*8 ZBSUI(KDLON)                    ! Intermediate variable
      REAL*8 ZCTS(KDLON,KFLEV)               ! Intermediate variable
      REAL*8 ZCNTRB(KDLON,KFLEV+1,KFLEV+1)   ! Intermediate variable
      SAVE ZFLUX, ZFLUC, ZBINT, ZBSUI, ZCTS, ZCNTRB
c
      INTEGER ilim, i, k, kpl1
C
      INTEGER lw0pas ! Every lw0pas steps, clear-sky is done
      PARAMETER (lw0pas=1)
      INTEGER lwpas  ! Every lwpas steps, cloudy-sky is done
      PARAMETER (lwpas=1)
c
      INTEGER itaplw0, itaplw
      LOGICAL appel1er
      SAVE appel1er, itaplw0, itaplw
      DATA appel1er /.TRUE./
      DATA itaplw0,itaplw /0,0/
C     ------------------------------------------------------------------
      IF (appel1er) THEN
         PRINT*, "LW clear-sky calling frequency: ", lw0pas
         PRINT*, "LW cloudy-sky calling frequency: ", lwpas
         PRINT*, "   In general, they should be 1"
         appel1er=.FALSE.
      ENDIF
C
      IF (MOD(itaplw0,lw0pas).EQ.0) THEN
      DO k = 1, KFLEV  ! convertir ozone de kg/kg en pa/pa
      DO i = 1, KDLON
c convertir ozone de kg/kg en pa (modif MPL 100505)
         ZOZ(i,k) = POZON(i,k)*PDP(i,k) * MD/RMO3
c        print *,'LW: ZOZ*10**6=',ZOZ(i,k)*1000000.
      ENDDO
      ENDDO
cIM ctes ds clesphys.h   CALL LWU(RCO2,RCH4, RN2O, RCFC11, RCFC12,
      CALL LWU(
     S         PAER,PDP,PPMB,PPSOL,ZOZ,PTAVE,PVIEW,PWV,ZABCU)
      CALL LWBV(ILIM,PDP,PDT0,PEMIS,PPMB,PTL,PTAVE,ZABCU,
     S          ZFLUC,ZBINT,ZBSUI,ZCTS,ZCNTRB)
      itaplw0 = 0
      ENDIF
      itaplw0 = itaplw0 + 1
C
      IF (MOD(itaplw,lwpas).EQ.0) THEN
      CALL LWC(ILIM,PCLDLD,PCLDLU,PEMIS,
     S         ZFLUC,ZBINT,ZBSUI,ZCTS,ZCNTRB,
     S         ZFLUX)
      itaplw = 0
      ENDIF
      itaplw = itaplw + 1
C
      DO k = 1, KFLEV
         kpl1 = k+1
         DO i = 1, KDLON
            PCOLR(i,k) = ZFLUX(i,1,kpl1)+ZFLUX(i,2,kpl1)
     .                 - ZFLUX(i,1,k)-   ZFLUX(i,2,k)
            PCOLR(i,k) = PCOLR(i,k) * RDAY*RG/RCPD / PDP(i,k)
            PCOLR0(i,k) = ZFLUC(i,1,kpl1)+ZFLUC(i,2,kpl1)
     .                 - ZFLUC(i,1,k)-   ZFLUC(i,2,k)
            PCOLR0(i,k) = PCOLR0(i,k) * RDAY*RG/RCPD / PDP(i,k)
         ENDDO
      ENDDO
      DO i = 1, KDLON
         PSOLLW(i) = -ZFLUX(i,1,1)-ZFLUX(i,2,1)
         PTOPLW(i) = ZFLUX(i,1,KFLEV+1) + ZFLUX(i,2,KFLEV+1)
c
         PSOLLW0(i) = -ZFLUC(i,1,1)-ZFLUC(i,2,1)
         PTOPLW0(i) = ZFLUC(i,1,KFLEV+1) + ZFLUC(i,2,KFLEV+1)
         psollwdown(i) = -ZFLUX(i,2,1)
c
cIM attention aux signes !; LWtop >0, LWdn < 0
         DO k = 1, KFLEV+1
           plwup(i,k) = ZFLUX(i,1,k)
           plwup0(i,k) = ZFLUC(i,1,k)
           plwdn(i,k) = ZFLUX(i,2,k)
           plwdn0(i,k) = ZFLUC(i,2,k)
         ENDDO
      ENDDO
C     ------------------------------------------------------------------
      RETURN
      END
