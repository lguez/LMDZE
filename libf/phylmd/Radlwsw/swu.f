c
cIM ctes ds clesphys.h   SUBROUTINE SWU (PSCT,RCO2,PCLDSW,PPMB,PPSOL,PRMU0,PFRAC,
      SUBROUTINE SWU (PSCT,PCLDSW,PPMB,PPSOL,PRMU0,PFRAC,
     S                PTAVE,PWV,PAKI,PCLD,PCLEAR,PDSIG,PFACT,
     S                PRMU,PSEC,PUD)
      use dimens_m
      use dimphy
      use clesphys
      use SUPHEC_M
      use raddim
      use radepsi
      use radopt
      IMPLICIT none
C
C* ARGUMENTS:
C
      REAL*8 PSCT
cIM ctes ds clesphys.h   REAL*8 RCO2
      REAL*8 PCLDSW(KDLON,KFLEV)
      REAL*8 PPMB(KDLON,KFLEV+1)
      REAL*8 PPSOL(KDLON)
      REAL*8 PRMU0(KDLON)
      REAL*8 PFRAC(KDLON)
      REAL*8 PTAVE(KDLON,KFLEV)
      REAL*8 PWV(KDLON,KFLEV)
C
      REAL*8 PAKI(KDLON,2)
      REAL*8 PCLD(KDLON,KFLEV)
      REAL*8 PCLEAR(KDLON)
      REAL*8 PDSIG(KDLON,KFLEV)
      REAL*8 PFACT(KDLON)
      REAL*8 PRMU(KDLON)
      REAL*8 PSEC(KDLON)
      REAL*8 PUD(KDLON,5,KFLEV+1)
C
C* LOCAL VARIABLES:
C
      INTEGER IIND(2)
      REAL*8 ZC1J(KDLON,KFLEV+1)
      REAL*8 ZCLEAR(KDLON)
      REAL*8 ZCLOUD(KDLON)
      REAL*8 ZN175(KDLON)
      REAL*8 ZN190(KDLON)
      REAL*8 ZO175(KDLON)
      REAL*8 ZO190(KDLON)
      REAL*8 ZSIGN(KDLON)
      REAL*8 ZR(KDLON,2) 
      REAL*8 ZSIGO(KDLON)
      REAL*8 ZUD(KDLON,2)
      REAL*8 ZRTH, ZRTU, ZWH2O, ZDSCO2, ZDSH2O, ZFPPW
      INTEGER jl, jk, jkp1, jkl, jklp1, ja
C
C* Prescribed Data:
c
      REAL*8 ZPDH2O,ZPDUMG
      SAVE ZPDH2O,ZPDUMG
      REAL*8 ZPRH2O,ZPRUMG
      SAVE ZPRH2O,ZPRUMG
      REAL*8 RTDH2O,RTDUMG
      SAVE RTDH2O,RTDUMG
      REAL*8 RTH2O ,RTUMG
      SAVE RTH2O ,RTUMG
      DATA ZPDH2O,ZPDUMG / 0.8   , 0.75 /
      DATA ZPRH2O,ZPRUMG / 30000., 30000. /
      DATA RTDH2O,RTDUMG /  0.40  , 0.375 /
      DATA RTH2O ,RTUMG  /  240.  , 240.  /
C     ------------------------------------------------------------------
C
C*         1.     COMPUTES AMOUNTS OF ABSORBERS
C                 -----------------------------
C
 100  CONTINUE
C
      IIND(1)=1
      IIND(2)=2
C      
C
C*         1.1    INITIALIZES QUANTITIES
C                 ----------------------
C
 110  CONTINUE
C
      DO 111 JL = 1, KDLON
      PUD(JL,1,KFLEV+1)=0.
      PUD(JL,2,KFLEV+1)=0.
      PUD(JL,3,KFLEV+1)=0.
      PUD(JL,4,KFLEV+1)=0.
      PUD(JL,5,KFLEV+1)=0.
      PFACT(JL)= PRMU0(JL) * PFRAC(JL) * PSCT
      PRMU(JL)=SQRT(1224.* PRMU0(JL) * PRMU0(JL) + 1.) / 35.
      PSEC(JL)=1./PRMU(JL)
      ZC1J(JL,KFLEV+1)=0.
 111  CONTINUE
C
C*          1.3    AMOUNTS OF ABSORBERS
C                  --------------------
C
 130  CONTINUE
C
      DO 131 JL= 1, KDLON
      ZUD(JL,1) = 0.
      ZUD(JL,2) = 0.
      ZO175(JL) = PPSOL(JL)** (ZPDUMG+1.)
      ZO190(JL) = PPSOL(JL)** (ZPDH2O+1.)
      ZSIGO(JL) = PPSOL(JL)
      ZCLEAR(JL)=1.
      ZCLOUD(JL)=0.
 131  CONTINUE
C
      DO 133 JK = 1 , KFLEV
      JKP1 = JK + 1
      JKL = KFLEV+1 - JK
      JKLP1 = JKL+1
      DO 132 JL = 1, KDLON
      ZRTH=(RTH2O/PTAVE(JL,JK))**RTDH2O
      ZRTU=(RTUMG/PTAVE(JL,JK))**RTDUMG
      ZWH2O = MAX (PWV(JL,JK) , ZEPSCQ )
      ZSIGN(JL) = 100. * PPMB(JL,JKP1)
      PDSIG(JL,JK) = (ZSIGO(JL) - ZSIGN(JL))/PPSOL(JL)
      ZN175(JL) = ZSIGN(JL) ** (ZPDUMG+1.)
      ZN190(JL) = ZSIGN(JL) ** (ZPDH2O+1.)
      ZDSCO2 = ZO175(JL) - ZN175(JL)
      ZDSH2O = ZO190(JL) - ZN190(JL)
      PUD(JL,1,JK) = 1./( 10.* RG * (ZPDH2O+1.) )/(ZPRH2O**ZPDH2O)
     .             * ZDSH2O * ZWH2O  * ZRTH
      PUD(JL,2,JK) = 1./( 10.* RG * (ZPDUMG+1.) )/(ZPRUMG**ZPDUMG)
     .             * ZDSCO2 * RCO2 * ZRTU
      ZFPPW=1.6078*ZWH2O/(1.+0.608*ZWH2O)
      PUD(JL,4,JK)=PUD(JL,1,JK)*ZFPPW
      PUD(JL,5,JK)=PUD(JL,1,JK)*(1.-ZFPPW)
      ZUD(JL,1) = ZUD(JL,1) + PUD(JL,1,JK)
      ZUD(JL,2) = ZUD(JL,2) + PUD(JL,2,JK)
      ZSIGO(JL) = ZSIGN(JL)
      ZO175(JL) = ZN175(JL)
      ZO190(JL) = ZN190(JL)
C      
      IF (NOVLP.EQ.1) THEN
         ZCLEAR(JL)=ZCLEAR(JL)
     S               *(1.-MAX(PCLDSW(JL,JKL),ZCLOUD(JL)))
     S               /(1.-MIN(ZCLOUD(JL),1.-ZEPSEC))
         ZC1J(JL,JKL)= 1.0 - ZCLEAR(JL)
         ZCLOUD(JL) = PCLDSW(JL,JKL)
      ELSE IF (NOVLP.EQ.2) THEN
         ZCLOUD(JL) = MAX(PCLDSW(JL,JKL),ZCLOUD(JL))
         ZC1J(JL,JKL) = ZCLOUD(JL)
      ELSE IF (NOVLP.EQ.3) THEN
         ZCLEAR(JL) = ZCLEAR(JL)*(1.-PCLDSW(JL,JKL))
         ZCLOUD(JL) = 1.0 - ZCLEAR(JL)
         ZC1J(JL,JKL) = ZCLOUD(JL)
      END IF
 132  CONTINUE
 133  CONTINUE
      DO 134 JL=1, KDLON
      PCLEAR(JL)=1.-ZC1J(JL,1)
 134  CONTINUE
      DO 136 JK=1,KFLEV
      DO 135 JL=1, KDLON
      IF (PCLEAR(JL).LT.1.) THEN
         PCLD(JL,JK)=PCLDSW(JL,JK)/(1.-PCLEAR(JL))
      ELSE
         PCLD(JL,JK)=0.
      END IF
 135  CONTINUE
 136  CONTINUE           
C      
C
C*         1.4    COMPUTES CLEAR-SKY GREY ABSORPTION COEFFICIENTS
C                 -----------------------------------------------
C
 140  CONTINUE
C
      DO 142 JA = 1,2
      DO 141 JL = 1, KDLON
      ZUD(JL,JA) = ZUD(JL,JA) * PSEC(JL)
 141  CONTINUE
 142  CONTINUE
C
      CALL SWTT1(2, 2, IIND, ZUD, ZR)
C
      DO 144 JA = 1,2
      DO 143 JL = 1, KDLON
      PAKI(JL,JA) = -LOG( ZR(JL,JA) ) / ZUD(JL,JA)
 143  CONTINUE
 144  CONTINUE
C
C
C     ------------------------------------------------------------------
C
      RETURN
      END
