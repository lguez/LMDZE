      SUBROUTINE SWR ( KNU
     S  , PALBD , PCG   , PCLD , PDSIG, POMEGA, PRAYL
     S  , PSEC  , PTAU
     S  , PCGAZ , PPIZAZ, PRAY1, PRAY2, PREFZ , PRJ  , PRK , PRMUE
     S  , PTAUAZ, PTRA1 , PTRA2 )
      use dimens_m
      use dimphy
      use raddim
      use radepsi
      use radopt
      IMPLICIT none
C
C     ------------------------------------------------------------------
C     PURPOSE.
C     --------
C           COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY IN CASE OF
C     CONTINUUM SCATTERING
C
C     METHOD.
C     -------
C
C          1. COMPUTES CONTINUUM FLUXES CORRESPONDING TO AEROSOL
C     OR/AND RAYLEIGH SCATTERING (NO MOLECULAR GAS ABSORPTION)
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
C     ------------------------------------------------------------------
C* ARGUMENTS:
C
      INTEGER KNU
      REAL*8 PALBD(KDLON,2)
      REAL*8 PCG(KDLON,2,KFLEV)
      REAL*8 PCLD(KDLON,KFLEV)
      REAL*8 PDSIG(KDLON,KFLEV)
      REAL*8 POMEGA(KDLON,2,KFLEV)
      REAL*8 PRAYL(KDLON)
      REAL*8 PSEC(KDLON)
      REAL*8 PTAU(KDLON,2,KFLEV)
C
      REAL*8 PRAY1(KDLON,KFLEV+1)
      REAL*8 PRAY2(KDLON,KFLEV+1)
      REAL*8 PREFZ(KDLON,2,KFLEV+1)
      REAL*8 PRJ(KDLON,6,KFLEV+1)
      REAL*8 PRK(KDLON,6,KFLEV+1)
      REAL*8 PRMUE(KDLON,KFLEV+1)
      REAL*8 PCGAZ(KDLON,KFLEV)
      REAL*8 PPIZAZ(KDLON,KFLEV)
      REAL*8 PTAUAZ(KDLON,KFLEV)
      REAL*8 PTRA1(KDLON,KFLEV+1)
      REAL*8 PTRA2(KDLON,KFLEV+1)
C
C* LOCAL VARIABLES:
C
      REAL*8 ZC1I(KDLON,KFLEV+1)
      REAL*8 ZCLEQ(KDLON,KFLEV)
      REAL*8 ZCLEAR(KDLON)
      REAL*8 ZCLOUD(KDLON)
      REAL*8 ZGG(KDLON)
      REAL*8 ZREF(KDLON)
      REAL*8 ZRE1(KDLON)
      REAL*8 ZRE2(KDLON)
      REAL*8 ZRMUZ(KDLON)
      REAL*8 ZRNEB(KDLON)
      REAL*8 ZR21(KDLON)
      REAL*8 ZR22(KDLON)
      REAL*8 ZR23(KDLON)
      REAL*8 ZSS1(KDLON)
      REAL*8 ZTO1(KDLON)
      REAL*8 ZTR(KDLON,2,KFLEV+1)
      REAL*8 ZTR1(KDLON)
      REAL*8 ZTR2(KDLON)
      REAL*8 ZW(KDLON)
C
      INTEGER jk, jl, ja, jkl, jklp1, jkm1, jaj
      REAL*8 ZFACOA, ZFACOC, ZCORAE, ZCORCD
      REAL*8 ZMUE, ZGAP, ZWW, ZTO, ZDEN, ZDEN1
      REAL*8 ZMU1, ZRE11, ZBMU0, ZBMU1
C
C     ------------------------------------------------------------------
C
C*         1.    INITIALIZATION
C                --------------
C
 100  CONTINUE
C
      DO 103 JK = 1 , KFLEV+1
      DO 102 JA = 1 , 6
      DO 101 JL = 1, KDLON
      PRJ(JL,JA,JK) = 0.
      PRK(JL,JA,JK) = 0.
 101  CONTINUE
 102  CONTINUE
 103  CONTINUE
C
C
C     ------------------------------------------------------------------
C
C*         2.    TOTAL EFFECTIVE CLOUDINESS ABOVE A GIVEN LEVEL
C                ----------------------------------------------
C
 200  CONTINUE
C
      DO 201 JL = 1, KDLON
      ZR23(JL) = 0.
      ZC1I(JL,KFLEV+1) = 0.
      ZCLEAR(JL) = 1.
      ZCLOUD(JL) = 0.
 201  CONTINUE
C
      JK = 1
      JKL = KFLEV+1 - JK
      JKLP1 = JKL + 1
      DO 202 JL = 1, KDLON
      ZFACOA = 1. - PPIZAZ(JL,JKL)*PCGAZ(JL,JKL)*PCGAZ(JL,JKL)
      ZFACOC = 1. - POMEGA(JL,KNU,JKL) * PCG(JL,KNU,JKL)
     S                                 * PCG(JL,KNU,JKL)
      ZCORAE = ZFACOA * PTAUAZ(JL,JKL) * PSEC(JL)
      ZCORCD = ZFACOC * PTAU(JL,KNU,JKL) * PSEC(JL)
      ZR21(JL) = EXP(-ZCORAE   )
      ZR22(JL) = EXP(-ZCORCD   )
      ZSS1(JL) = PCLD(JL,JKL)*(1.0-ZR21(JL)*ZR22(JL))
     S               + (1.0-PCLD(JL,JKL))*(1.0-ZR21(JL))
      ZCLEQ(JL,JKL) = ZSS1(JL)
C
      IF (NOVLP.EQ.1) THEN
c* maximum-random
         ZCLEAR(JL) = ZCLEAR(JL)
     S                  *(1.0-MAX(ZSS1(JL),ZCLOUD(JL)))
     S                  /(1.0-MIN(ZCLOUD(JL),1.-ZEPSEC))
         ZC1I(JL,JKL) = 1.0 - ZCLEAR(JL)
         ZCLOUD(JL) = ZSS1(JL)
      ELSE IF (NOVLP.EQ.2) THEN
C* maximum
         ZCLOUD(JL) = MAX( ZSS1(JL) , ZCLOUD(JL) )
         ZC1I(JL,JKL) = ZCLOUD(JL)
      ELSE IF (NOVLP.EQ.3) THEN
c* random
         ZCLEAR(JL) = ZCLEAR(JL)*(1.0 - ZSS1(JL))
         ZCLOUD(JL) = 1.0 - ZCLEAR(JL)
         ZC1I(JL,JKL) = ZCLOUD(JL)
      END IF
 202  CONTINUE
C
      DO 205 JK = 2 , KFLEV
      JKL = KFLEV+1 - JK
      JKLP1 = JKL + 1
      DO 204 JL = 1, KDLON
      ZFACOA = 1. - PPIZAZ(JL,JKL)*PCGAZ(JL,JKL)*PCGAZ(JL,JKL)
      ZFACOC = 1. - POMEGA(JL,KNU,JKL) * PCG(JL,KNU,JKL)
     S                                 * PCG(JL,KNU,JKL)
      ZCORAE = ZFACOA * PTAUAZ(JL,JKL) * PSEC(JL)
      ZCORCD = ZFACOC * PTAU(JL,KNU,JKL) * PSEC(JL)
      ZR21(JL) = EXP(-ZCORAE   )
      ZR22(JL) = EXP(-ZCORCD   )
      ZSS1(JL) = PCLD(JL,JKL)*(1.0-ZR21(JL)*ZR22(JL))
     S               + (1.0-PCLD(JL,JKL))*(1.0-ZR21(JL))
      ZCLEQ(JL,JKL) = ZSS1(JL)
c     
      IF (NOVLP.EQ.1) THEN
c* maximum-random
         ZCLEAR(JL) = ZCLEAR(JL)
     S                  *(1.0-MAX(ZSS1(JL),ZCLOUD(JL)))
     S                  /(1.0-MIN(ZCLOUD(JL),1.-ZEPSEC))
         ZC1I(JL,JKL) = 1.0 - ZCLEAR(JL)
         ZCLOUD(JL) = ZSS1(JL)
      ELSE IF (NOVLP.EQ.2) THEN
C* maximum
         ZCLOUD(JL) = MAX( ZSS1(JL) , ZCLOUD(JL) )
         ZC1I(JL,JKL) = ZCLOUD(JL)
      ELSE IF (NOVLP.EQ.3) THEN
c* random
         ZCLEAR(JL) = ZCLEAR(JL)*(1.0 - ZSS1(JL))
         ZCLOUD(JL) = 1.0 - ZCLEAR(JL)
         ZC1I(JL,JKL) = ZCLOUD(JL)
      END IF
 204  CONTINUE
 205  CONTINUE
C
C     ------------------------------------------------------------------
C
C*         3.    REFLECTIVITY/TRANSMISSIVITY FOR PURE SCATTERING
C                -----------------------------------------------
C
 300  CONTINUE
C
      DO 301 JL = 1, KDLON
      PRAY1(JL,KFLEV+1) = 0.
      PRAY2(JL,KFLEV+1) = 0.
      PREFZ(JL,2,1) = PALBD(JL,KNU)
      PREFZ(JL,1,1) = PALBD(JL,KNU)
      PTRA1(JL,KFLEV+1) = 1.
      PTRA2(JL,KFLEV+1) = 1.
 301  CONTINUE
C
      DO 346 JK = 2 , KFLEV+1
      JKM1 = JK-1
      DO 342 JL = 1, KDLON
      ZRNEB(JL)= PCLD(JL,JKM1)
      ZRE1(JL)=0.
      ZTR1(JL)=0.
      ZRE2(JL)=0.
      ZTR2(JL)=0.
C
C
C     ------------------------------------------------------------------
C
C*         3.1  EQUIVALENT ZENITH ANGLE
C               -----------------------
C
 310  CONTINUE
C
      ZMUE = (1.-ZC1I(JL,JK)) * PSEC(JL)
     S            + ZC1I(JL,JK) * 1.66
      PRMUE(JL,JK) = 1./ZMUE
C
C
C     ------------------------------------------------------------------
C
C*         3.2  REFLECT./TRANSMISSIVITY DUE TO RAYLEIGH AND AEROSOLS
C               ----------------------------------------------------
C
 320  CONTINUE
C
      ZGAP = PCGAZ(JL,JKM1)
      ZBMU0 = 0.5 - 0.75 * ZGAP / ZMUE
      ZWW = PPIZAZ(JL,JKM1)
      ZTO = PTAUAZ(JL,JKM1)
      ZDEN = 1. + (1. - ZWW + ZBMU0 * ZWW) * ZTO * ZMUE
     S       + (1-ZWW) * (1. - ZWW +2.*ZBMU0*ZWW)*ZTO*ZTO*ZMUE*ZMUE
      PRAY1(JL,JKM1) = ZBMU0 * ZWW * ZTO * ZMUE / ZDEN
      PTRA1(JL,JKM1) = 1. / ZDEN
c      PRINT *,' LOOP 342 ** 3 ** JL=',JL,PRAY1(JL,JKM1),PTRA1(JL,JKM1)
C
      ZMU1 = 0.5
      ZBMU1 = 0.5 - 0.75 * ZGAP * ZMU1
      ZDEN1= 1. + (1. - ZWW + ZBMU1 * ZWW) * ZTO / ZMU1
     S       + (1-ZWW) * (1. - ZWW +2.*ZBMU1*ZWW)*ZTO*ZTO/ZMU1/ZMU1
      PRAY2(JL,JKM1) = ZBMU1 * ZWW * ZTO / ZMU1 / ZDEN1
      PTRA2(JL,JKM1) = 1. / ZDEN1
C
C
C     ------------------------------------------------------------------
C
C*         3.3  EFFECT OF CLOUD LAYER
C               ---------------------
C
 330  CONTINUE
C
      ZW(JL) = POMEGA(JL,KNU,JKM1)
      ZTO1(JL) = PTAU(JL,KNU,JKM1)/ZW(JL)
     S         + PTAUAZ(JL,JKM1)/PPIZAZ(JL,JKM1)
      ZR21(JL) = PTAU(JL,KNU,JKM1) + PTAUAZ(JL,JKM1)
      ZR22(JL) = PTAU(JL,KNU,JKM1) / ZR21(JL)
      ZGG(JL) = ZR22(JL) * PCG(JL,KNU,JKM1)
     S              + (1. - ZR22(JL)) * PCGAZ(JL,JKM1)
C Modif PhD - JJM 19/03/96 pour erreurs arrondis
C machine
C PHD PROTECTION ZW(JL) = ZR21(JL) / ZTO1(JL)
      IF (ZW(JL).EQ.1. .AND. PPIZAZ(JL,JKM1).EQ.1.) THEN
         ZW(JL)=1.
      ELSE
         ZW(JL) = ZR21(JL) / ZTO1(JL)
      END IF
      ZREF(JL) = PREFZ(JL,1,JKM1)
      ZRMUZ(JL) = PRMUE(JL,JK)
 342  CONTINUE
C
      CALL SWDE(ZGG  , ZREF  , ZRMUZ , ZTO1 , ZW,
     S          ZRE1 , ZRE2  , ZTR1  , ZTR2)
C
      DO 345 JL = 1, KDLON
C
      PREFZ(JL,1,JK) = (1.-ZRNEB(JL)) * (PRAY1(JL,JKM1)
     S               + PREFZ(JL,1,JKM1) * PTRA1(JL,JKM1)
     S               * PTRA2(JL,JKM1)
     S               / (1.-PRAY2(JL,JKM1)*PREFZ(JL,1,JKM1)))
     S               + ZRNEB(JL) * ZRE2(JL)
C
      ZTR(JL,1,JKM1) = ZRNEB(JL) * ZTR2(JL) + (PTRA1(JL,JKM1)
     S               / (1.-PRAY2(JL,JKM1)*PREFZ(JL,1,JKM1)))
     S               * (1.-ZRNEB(JL))
C
      PREFZ(JL,2,JK) = (1.-ZRNEB(JL)) * (PRAY1(JL,JKM1)
     S               + PREFZ(JL,2,JKM1) * PTRA1(JL,JKM1)
     S               * PTRA2(JL,JKM1) )
     S               + ZRNEB(JL) * ZRE1(JL)
C
      ZTR(JL,2,JKM1) = ZRNEB(JL) * ZTR1(JL)
     S               + PTRA1(JL,JKM1) * (1.-ZRNEB(JL))
C
 345  CONTINUE
 346  CONTINUE
      DO 347 JL = 1, KDLON
      ZMUE = (1.-ZC1I(JL,1))*PSEC(JL)+ZC1I(JL,1)*1.66
      PRMUE(JL,1)=1./ZMUE
 347  CONTINUE
C
C
C     ------------------------------------------------------------------
C
C*         3.5    REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL
C                 -------------------------------------------------
C
 350  CONTINUE
C
      IF (KNU.EQ.1) THEN
      JAJ = 2
      DO 351 JL = 1, KDLON
      PRJ(JL,JAJ,KFLEV+1) = 1.
      PRK(JL,JAJ,KFLEV+1) = PREFZ(JL, 1,KFLEV+1)
 351  CONTINUE
C
      DO 353 JK = 1 , KFLEV
      JKL = KFLEV+1 - JK
      JKLP1 = JKL + 1
      DO 352 JL = 1, KDLON
      ZRE11= PRJ(JL,JAJ,JKLP1) * ZTR(JL,  1,JKL)
      PRJ(JL,JAJ,JKL) = ZRE11
      PRK(JL,JAJ,JKL) = ZRE11 * PREFZ(JL,  1,JKL)
 352  CONTINUE
 353  CONTINUE
 354  CONTINUE
C
      ELSE
C
      DO 358 JAJ = 1 , 2
      DO 355 JL = 1, KDLON
      PRJ(JL,JAJ,KFLEV+1) = 1.
      PRK(JL,JAJ,KFLEV+1) = PREFZ(JL,JAJ,KFLEV+1)
 355  CONTINUE
C
      DO 357 JK = 1 , KFLEV
      JKL = KFLEV+1 - JK
      JKLP1 = JKL + 1
      DO 356 JL = 1, KDLON
      ZRE11= PRJ(JL,JAJ,JKLP1) * ZTR(JL,JAJ,JKL)
      PRJ(JL,JAJ,JKL) = ZRE11
      PRK(JL,JAJ,JKL) = ZRE11 * PREFZ(JL,JAJ,JKL)
 356  CONTINUE
 357  CONTINUE
 358  CONTINUE
C
      END IF
C
C     ------------------------------------------------------------------
C
      RETURN
      END
