      SUBROUTINE SWCLR  ( KNU
     S  , PAER  , flag_aer, tauae, pizae, cgae
     S  , PALBP , PDSIG , PRAYL , PSEC
     S  , PCGAZ , PPIZAZ, PRAY1 , PRAY2 , PREFZ , PRJ  
     S  , PRK   , PRMU0 , PTAUAZ, PTRA1 , PTRA2                   )
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
C     CLEAR-SKY COLUMN
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
C        ORIGINAL : 94-11-15
C     ------------------------------------------------------------------
C* ARGUMENTS:
C
      INTEGER KNU
c-OB
      real*8 flag_aer
      real*8 tauae(kdlon,kflev,2)
      real*8 pizae(kdlon,kflev,2)
      real*8 cgae(kdlon,kflev,2)
      REAL*8 PAER(KDLON,KFLEV,5)
      REAL*8 PALBP(KDLON,2)
      REAL*8 PDSIG(KDLON,KFLEV)
      REAL*8 PRAYL(KDLON)
      REAL*8 PSEC(KDLON)
C
      REAL*8 PCGAZ(KDLON,KFLEV)     
      REAL*8 PPIZAZ(KDLON,KFLEV)
      REAL*8 PRAY1(KDLON,KFLEV+1)
      REAL*8 PRAY2(KDLON,KFLEV+1)
      REAL*8 PREFZ(KDLON,2,KFLEV+1)
      REAL*8 PRJ(KDLON,6,KFLEV+1)
      REAL*8 PRK(KDLON,6,KFLEV+1)
      REAL*8 PRMU0(KDLON,KFLEV+1)
      REAL*8 PTAUAZ(KDLON,KFLEV)
      REAL*8 PTRA1(KDLON,KFLEV+1)
      REAL*8 PTRA2(KDLON,KFLEV+1)
C
C* LOCAL VARIABLES:
C
      REAL*8 ZC0I(KDLON,KFLEV+1)       
      REAL*8 ZCLE0(KDLON,KFLEV)
      REAL*8 ZCLEAR(KDLON)
      REAL*8 ZR21(KDLON)
      REAL*8 ZR23(KDLON)
      REAL*8 ZSS0(KDLON)
      REAL*8 ZSCAT(KDLON)
      REAL*8 ZTR(KDLON,2,KFLEV+1)
C
      INTEGER jl, jk, ja, jae, jkl, jklp1, jaj, jkm1, in
      REAL*8 ZTRAY, ZGAR, ZRATIO, ZFF, ZFACOA, ZCORAE
      REAL*8 ZMUE, ZGAP, ZWW, ZTO, ZDEN, ZMU1, ZDEN1
      REAL*8 ZBMU0, ZBMU1, ZRE11
C
C* Prescribed Data for Aerosols:
C
      REAL*8 TAUA(2,5), RPIZA(2,5), RCGA(2,5)
      SAVE TAUA, RPIZA, RCGA
      DATA ((TAUA(IN,JA),JA=1,5),IN=1,2) /
     S .730719, .912819, .725059, .745405, .682188 ,
     S .730719, .912819, .725059, .745405, .682188 /
      DATA ((RPIZA(IN,JA),JA=1,5),IN=1,2) /
     S .872212, .982545, .623143, .944887, .997975 ,
     S .872212, .982545, .623143, .944887, .997975 /
      DATA ((RCGA (IN,JA),JA=1,5),IN=1,2) /
     S .647596, .739002, .580845, .662657, .624246 ,
     S .647596, .739002, .580845, .662657, .624246 /
C     ------------------------------------------------------------------
C
C*         1.    OPTICAL PARAMETERS FOR AEROSOLS AND RAYLEIGH
C                --------------------------------------------
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
      DO 108 JK = 1 , KFLEV
c-OB
c      DO 104 JL = 1, KDLON
c      PCGAZ(JL,JK) = 0.
c      PPIZAZ(JL,JK) =  0.
c      PTAUAZ(JL,JK) = 0.
c 104  CONTINUE
c-OB
c      DO 106 JAE=1,5
c      DO 105 JL = 1, KDLON
c      PTAUAZ(JL,JK)=PTAUAZ(JL,JK)
c     S        +PAER(JL,JK,JAE)*TAUA(KNU,JAE)
c      PPIZAZ(JL,JK)=PPIZAZ(JL,JK)+PAER(JL,JK,JAE)
c     S        * TAUA(KNU,JAE)*RPIZA(KNU,JAE)
c      PCGAZ(JL,JK) =  PCGAZ(JL,JK) +PAER(JL,JK,JAE)
c     S        * TAUA(KNU,JAE)*RPIZA(KNU,JAE)*RCGA(KNU,JAE)
c 105  CONTINUE
c 106  CONTINUE
c-OB
      DO 105 JL = 1, KDLON
      PTAUAZ(JL,JK)=flag_aer * tauae(JL,JK,KNU)
      PPIZAZ(JL,JK)=flag_aer * pizae(JL,JK,KNU)
      PCGAZ (JL,JK)=flag_aer * cgae(JL,JK,KNU)
 105  CONTINUE
C
      IF (flag_aer.GT.0) THEN
c-OB
      DO 107 JL = 1, KDLON
c         PCGAZ(JL,JK)=PCGAZ(JL,JK)/PPIZAZ(JL,JK)
c         PPIZAZ(JL,JK)=PPIZAZ(JL,JK)/PTAUAZ(JL,JK)
         ZTRAY = PRAYL(JL) * PDSIG(JL,JK)
         ZRATIO = ZTRAY / (ZTRAY + PTAUAZ(JL,JK))
         ZGAR = PCGAZ(JL,JK)
         ZFF = ZGAR * ZGAR
         PTAUAZ(JL,JK)=ZTRAY+PTAUAZ(JL,JK)*(1.-PPIZAZ(JL,JK)*ZFF)
         PCGAZ(JL,JK) = ZGAR * (1. - ZRATIO) / (1. + ZGAR)
         PPIZAZ(JL,JK) =ZRATIO+(1.-ZRATIO)*PPIZAZ(JL,JK)*(1.-ZFF)
     S                       / (1. - PPIZAZ(JL,JK) * ZFF)
 107  CONTINUE
      ELSE
      DO JL = 1, KDLON
         ZTRAY = PRAYL(JL) * PDSIG(JL,JK)
         PTAUAZ(JL,JK) = ZTRAY
         PCGAZ(JL,JK) = 0.
         PPIZAZ(JL,JK) = 1.-REPSCT
      END DO
      END IF   ! check flag_aer
c     107  CONTINUE
c      PRINT 9107,JK,((PAER(JL,JK,JAE),JAE=1,5)
c     $ ,PTAUAZ(JL,JK),PPIZAZ(JL,JK),PCGAZ(JL,JK),JL=1,KDLON)
c 9107 FORMAT(1X,'SWCLR_107',I3,8E12.5)
C
 108  CONTINUE
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
      ZC0I(JL,KFLEV+1) = 0.
      ZCLEAR(JL) = 1.
      ZSCAT(JL) = 0.
 201  CONTINUE
C
      JK = 1
      JKL = KFLEV+1 - JK
      JKLP1 = JKL + 1
      DO 202 JL = 1, KDLON
      ZFACOA = 1. - PPIZAZ(JL,JKL)*PCGAZ(JL,JKL)*PCGAZ(JL,JKL)
      ZCORAE = ZFACOA * PTAUAZ(JL,JKL) * PSEC(JL)
      ZR21(JL) = EXP(-ZCORAE   )
      ZSS0(JL) = 1.-ZR21(JL)
      ZCLE0(JL,JKL) = ZSS0(JL)
C
      IF (NOVLP.EQ.1) THEN
c* maximum-random
         ZCLEAR(JL) = ZCLEAR(JL)
     S                  *(1.0-MAX(ZSS0(JL),ZSCAT(JL)))
     S                  /(1.0-MIN(ZSCAT(JL),1.-ZEPSEC))
         ZC0I(JL,JKL) = 1.0 - ZCLEAR(JL)
         ZSCAT(JL) = ZSS0(JL)
      ELSE IF (NOVLP.EQ.2) THEN
C* maximum
         ZSCAT(JL) = MAX( ZSS0(JL) , ZSCAT(JL) )
         ZC0I(JL,JKL) = ZSCAT(JL)
      ELSE IF (NOVLP.EQ.3) THEN
c* random
         ZCLEAR(JL)=ZCLEAR(JL)*(1.0-ZSS0(JL))
         ZSCAT(JL) = 1.0 - ZCLEAR(JL)
         ZC0I(JL,JKL) = ZSCAT(JL)
      END IF
 202  CONTINUE
C
      DO 205 JK = 2 , KFLEV
      JKL = KFLEV+1 - JK
      JKLP1 = JKL + 1
      DO 204 JL = 1, KDLON
      ZFACOA = 1. - PPIZAZ(JL,JKL)*PCGAZ(JL,JKL)*PCGAZ(JL,JKL)
      ZCORAE = ZFACOA * PTAUAZ(JL,JKL) * PSEC(JL)
      ZR21(JL) = EXP(-ZCORAE   )
      ZSS0(JL) = 1.-ZR21(JL)
      ZCLE0(JL,JKL) = ZSS0(JL)
c     
      IF (NOVLP.EQ.1) THEN
c* maximum-random
         ZCLEAR(JL) = ZCLEAR(JL)
     S                  *(1.0-MAX(ZSS0(JL),ZSCAT(JL)))
     S                  /(1.0-MIN(ZSCAT(JL),1.-ZEPSEC))
         ZC0I(JL,JKL) = 1.0 - ZCLEAR(JL)
         ZSCAT(JL) = ZSS0(JL)
      ELSE IF (NOVLP.EQ.2) THEN
C* maximum
         ZSCAT(JL) = MAX( ZSS0(JL) , ZSCAT(JL) )
         ZC0I(JL,JKL) = ZSCAT(JL)
      ELSE IF (NOVLP.EQ.3) THEN
c* random
         ZCLEAR(JL)=ZCLEAR(JL)*(1.0-ZSS0(JL))
         ZSCAT(JL) = 1.0 - ZCLEAR(JL)
         ZC0I(JL,JKL) = ZSCAT(JL)
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
      PREFZ(JL,2,1) = PALBP(JL,KNU)
      PREFZ(JL,1,1) = PALBP(JL,KNU)
      PTRA1(JL,KFLEV+1) = 1.
      PTRA2(JL,KFLEV+1) = 1.
 301  CONTINUE
C
      DO 346 JK = 2 , KFLEV+1
      JKM1 = JK-1
      DO 342 JL = 1, KDLON
C
C
C     ------------------------------------------------------------------
C
C*         3.1  EQUIVALENT ZENITH ANGLE
C               -----------------------
C
 310  CONTINUE
C
      ZMUE = (1.-ZC0I(JL,JK)) * PSEC(JL)
     S            + ZC0I(JL,JK) * 1.66
      PRMU0(JL,JK) = 1./ZMUE
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
C
      ZMU1 = 0.5
      ZBMU1 = 0.5 - 0.75 * ZGAP * ZMU1
      ZDEN1= 1. + (1. - ZWW + ZBMU1 * ZWW) * ZTO / ZMU1
     S       + (1-ZWW) * (1. - ZWW +2.*ZBMU1*ZWW)*ZTO*ZTO/ZMU1/ZMU1
      PRAY2(JL,JKM1) = ZBMU1 * ZWW * ZTO / ZMU1 / ZDEN1
      PTRA2(JL,JKM1) = 1. / ZDEN1
C
C
C
      PREFZ(JL,1,JK) = (PRAY1(JL,JKM1)
     S               + PREFZ(JL,1,JKM1) * PTRA1(JL,JKM1)
     S               * PTRA2(JL,JKM1)
     S               / (1.-PRAY2(JL,JKM1)*PREFZ(JL,1,JKM1)))
C
      ZTR(JL,1,JKM1) = (PTRA1(JL,JKM1)
     S               / (1.-PRAY2(JL,JKM1)*PREFZ(JL,1,JKM1)))
C
      PREFZ(JL,2,JK) = (PRAY1(JL,JKM1)
     S               + PREFZ(JL,2,JKM1) * PTRA1(JL,JKM1)
     S               * PTRA2(JL,JKM1) )
C
      ZTR(JL,2,JKM1) = PTRA1(JL,JKM1) 
C
 342  CONTINUE
 346  CONTINUE
      DO 347 JL = 1, KDLON
      ZMUE = (1.-ZC0I(JL,1))*PSEC(JL)+ZC0I(JL,1)*1.66
      PRMU0(JL,1)=1./ZMUE
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
