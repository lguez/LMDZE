      SUBROUTINE SW2S ( KNU
     S  ,  PAER  , flag_aer, tauae, pizae, cgae
     S  ,  PAKI, PALBD, PALBP, PCG   , PCLD, PCLEAR, PCLDSW
     S  ,  PDSIG ,POMEGA,POZ , PRMU , PSEC  , PTAU
     S  ,  PUD   ,PWV , PQS
     S  ,  PFDOWN,PFUP                                            )
      use dimens_m
      use dimphy
      use raddim
      use radepsi
      IMPLICIT none
C
C     ------------------------------------------------------------------
C     PURPOSE.
C     --------
C
C          THIS ROUTINE COMPUTES THE SHORTWAVE RADIATION FLUXES IN THE
C     SECOND SPECTRAL INTERVAL FOLLOWING FOUQUART AND BONNEL (1980).
C
C     METHOD.
C     -------
C
C          1. COMPUTES REFLECTIVITY/TRANSMISSIVITY CORRESPONDING TO
C     CONTINUUM SCATTERING
C          2. COMPUTES REFLECTIVITY/TRANSMISSIVITY CORRESPONDING FOR
C     A GREY MOLECULAR ABSORPTION
C          3. LAPLACE TRANSFORM ON THE PREVIOUS TO GET EFFECTIVE AMOUNTS
C     OF ABSORBERS
C          4. APPLY H2O AND U.M.G. TRANSMISSION FUNCTIONS
C          5. MULTIPLY BY OZONE TRANSMISSION FUNCTION
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
C        94-11-15   J.-J. MORCRETTE    DIRECT/DIFFUSE ALBEDO
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
      REAL*8 PAKI(KDLON,2)
      REAL*8 PALBD(KDLON,2)
      REAL*8 PALBP(KDLON,2)
      REAL*8 PCG(KDLON,2,KFLEV)
      REAL*8 PCLD(KDLON,KFLEV)
      REAL*8 PCLDSW(KDLON,KFLEV)
      REAL*8 PCLEAR(KDLON)
      REAL*8 PDSIG(KDLON,KFLEV)
      REAL*8 POMEGA(KDLON,2,KFLEV)
      REAL*8 POZ(KDLON,KFLEV)
      REAL*8 PQS(KDLON,KFLEV)
      REAL*8 PRMU(KDLON)
      REAL*8 PSEC(KDLON)
      REAL*8 PTAU(KDLON,2,KFLEV)
      REAL*8 PUD(KDLON,5,KFLEV+1)
      REAL*8 PWV(KDLON,KFLEV)
C
      REAL*8 PFDOWN(KDLON,KFLEV+1)
      REAL*8 PFUP(KDLON,KFLEV+1)
C
C* LOCAL VARIABLES:
C
      INTEGER IIND2(2), IIND3(3)
      REAL*8 ZCGAZ(KDLON,KFLEV)
      REAL*8 ZFD(KDLON,KFLEV+1)
      REAL*8 ZFU(KDLON,KFLEV+1) 
      REAL*8 ZG(KDLON)
      REAL*8 ZGG(KDLON)
      REAL*8 ZPIZAZ(KDLON,KFLEV)
      REAL*8 ZRAYL(KDLON)
      REAL*8 ZRAY1(KDLON,KFLEV+1)
      REAL*8 ZRAY2(KDLON,KFLEV+1)
      REAL*8 ZREF(KDLON)
      REAL*8 ZREFZ(KDLON,2,KFLEV+1)
      REAL*8 ZRE1(KDLON)
      REAL*8 ZRE2(KDLON)
      REAL*8 ZRJ(KDLON,6,KFLEV+1)
      REAL*8 ZRJ0(KDLON,6,KFLEV+1)
      REAL*8 ZRK(KDLON,6,KFLEV+1)
      REAL*8 ZRK0(KDLON,6,KFLEV+1)
      REAL*8 ZRL(KDLON,8)
      REAL*8 ZRMUE(KDLON,KFLEV+1)
      REAL*8 ZRMU0(KDLON,KFLEV+1)
      REAL*8 ZRMUZ(KDLON)
      REAL*8 ZRNEB(KDLON)
      REAL*8 ZRUEF(KDLON,8)
      REAL*8 ZR1(KDLON) 
      REAL*8 ZR2(KDLON,2)
      REAL*8 ZR3(KDLON,3)
      REAL*8 ZR4(KDLON)
      REAL*8 ZR21(KDLON)
      REAL*8 ZR22(KDLON)
      REAL*8 ZS(KDLON)
      REAL*8 ZTAUAZ(KDLON,KFLEV)
      REAL*8 ZTO1(KDLON)
      REAL*8 ZTR(KDLON,2,KFLEV+1)
      REAL*8 ZTRA1(KDLON,KFLEV+1)
      REAL*8 ZTRA2(KDLON,KFLEV+1)
      REAL*8 ZTR1(KDLON)
      REAL*8 ZTR2(KDLON)
      REAL*8 ZW(KDLON)   
      REAL*8 ZW1(KDLON)
      REAL*8 ZW2(KDLON,2)
      REAL*8 ZW3(KDLON,3)
      REAL*8 ZW4(KDLON)
      REAL*8 ZW5(KDLON)
C
      INTEGER jl, jk, k, jaj, ikm1, ikl, jn, jabs, jkm1
      INTEGER jref, jkl, jklp1, jajp, jkki, jkkp4, jn2j, iabs
      REAL*8 ZRMUM1, ZWH2O, ZCNEB, ZAA, ZBB, ZRKI, ZRE11
C
C* Prescribed Data:
C
      REAL*8 RSUN(2)
      SAVE RSUN
      REAL*8 RRAY(2,6)
      SAVE RRAY
      DATA RSUN(1) / 0.441676 /
      DATA RSUN(2) / 0.558324 /
      DATA (RRAY(1,K),K=1,6) /
     S .428937E-01, .890743E+00,-.288555E+01,
     S .522744E+01,-.469173E+01, .161645E+01/
      DATA (RRAY(2,K),K=1,6) /
     S .697200E-02, .173297E-01,-.850903E-01,
     S .248261E+00,-.302031E+00, .129662E+00/
C
C     ------------------------------------------------------------------
C
C*         1.     SECOND SPECTRAL INTERVAL (0.68-4.00 MICRON)
C                 -------------------------------------------
C
 100  CONTINUE
C
C
C*         1.1    OPTICAL THICKNESS FOR RAYLEIGH SCATTERING
C                 -----------------------------------------
C
 110  CONTINUE
C
      DO 111 JL = 1, KDLON
      ZRMUM1 = 1. - PRMU(JL)
      ZRAYL(JL) =  RRAY(KNU,1) + ZRMUM1   * (RRAY(KNU,2) + ZRMUM1
     S          * (RRAY(KNU,3) + ZRMUM1   * (RRAY(KNU,4) + ZRMUM1
     S          * (RRAY(KNU,5) + ZRMUM1   *  RRAY(KNU,6)     ))))
 111  CONTINUE
C
C
C     ------------------------------------------------------------------
C
C*         2.    CONTINUUM SCATTERING CALCULATIONS
C                ---------------------------------
C
 200  CONTINUE
C
C*         2.1   CLEAR-SKY FRACTION OF THE COLUMN
C                --------------------------------
C  
 210  CONTINUE
C
      CALL SWCLR ( KNU
     S  , PAER   , flag_aer, tauae, pizae, cgae
     S  , PALBP  , PDSIG , ZRAYL, PSEC 
     S  , ZCGAZ  , ZPIZAZ, ZRAY1 , ZRAY2, ZREFZ, ZRJ0
     S  , ZRK0   , ZRMU0 , ZTAUAZ, ZTRA1, ZTRA2)
C
C
C*         2.2   CLOUDY FRACTION OF THE COLUMN
C                -----------------------------
C
 220  CONTINUE
C
      CALL SWR ( KNU
     S  , PALBD , PCG   , PCLD , PDSIG, POMEGA, ZRAYL
     S  , PSEC  , PTAU
     S  , ZCGAZ , ZPIZAZ, ZRAY1, ZRAY2, ZREFZ , ZRJ  , ZRK, ZRMUE
     S  , ZTAUAZ, ZTRA1 , ZTRA2)
C
C
C     ------------------------------------------------------------------
C
C*         3.    SCATTERING CALCULATIONS WITH GREY MOLECULAR ABSORPTION
C                ------------------------------------------------------
C
 300  CONTINUE
C
      JN = 2
C
      DO 361 JABS=1,2
C
C
C*         3.1  SURFACE CONDITIONS
C               ------------------
C
 310  CONTINUE
C
      DO 311 JL = 1, KDLON
      ZREFZ(JL,2,1) = PALBD(JL,KNU)
      ZREFZ(JL,1,1) = PALBD(JL,KNU)
 311  CONTINUE
C
C
C*         3.2  INTRODUCING CLOUD EFFECTS
C               -------------------------
C
 320  CONTINUE
C
      DO 324 JK = 2 , KFLEV+1
      JKM1 = JK - 1
      IKL=KFLEV+1-JKM1
      DO 322 JL = 1, KDLON
      ZRNEB(JL) = PCLD(JL,JKM1)
      IF (JABS.EQ.1 .AND. ZRNEB(JL).GT.2.*ZEELOG) THEN
         ZWH2O=MAX(PWV(JL,JKM1),ZEELOG)
         ZCNEB=MAX(ZEELOG,MIN(ZRNEB(JL),1.-ZEELOG))
         ZBB=PUD(JL,JABS,JKM1)*PQS(JL,JKM1)/ZWH2O
         ZAA=MAX((PUD(JL,JABS,JKM1)-ZCNEB*ZBB)/(1.-ZCNEB),ZEELOG)
      ELSE
         ZAA=PUD(JL,JABS,JKM1)
         ZBB=ZAA
      END IF
      ZRKI = PAKI(JL,JABS)
      ZS(JL) = EXP(-ZRKI * ZAA * 1.66)
      ZG(JL) = EXP(-ZRKI * ZAA / ZRMUE(JL,JK))
      ZTR1(JL) = 0.
      ZRE1(JL) = 0.
      ZTR2(JL) = 0.
      ZRE2(JL) = 0.
C
      ZW(JL)= POMEGA(JL,KNU,JKM1)
      ZTO1(JL) = PTAU(JL,KNU,JKM1) / ZW(JL)
     S               + ZTAUAZ(JL,JKM1) / ZPIZAZ(JL,JKM1)
     S               + ZBB * ZRKI

      ZR21(JL) = PTAU(JL,KNU,JKM1) + ZTAUAZ(JL,JKM1)
      ZR22(JL) = PTAU(JL,KNU,JKM1) / ZR21(JL)
      ZGG(JL) = ZR22(JL) * PCG(JL,KNU,JKM1)
     S              + (1. - ZR22(JL)) * ZCGAZ(JL,JKM1)
      ZW(JL) = ZR21(JL) / ZTO1(JL)
      ZREF(JL) = ZREFZ(JL,1,JKM1)
      ZRMUZ(JL) = ZRMUE(JL,JK)
 322  CONTINUE
C
      CALL SWDE(ZGG, ZREF, ZRMUZ, ZTO1, ZW,
     S          ZRE1, ZRE2, ZTR1, ZTR2)
C
      DO 323 JL = 1, KDLON
C
      ZREFZ(JL,2,JK) = (1.-ZRNEB(JL)) * (ZRAY1(JL,JKM1)
     S               + ZREFZ(JL,2,JKM1) * ZTRA1(JL,JKM1)
     S               * ZTRA2(JL,JKM1) ) * ZG(JL) * ZS(JL)
     S               + ZRNEB(JL) * ZRE1(JL)
C
      ZTR(JL,2,JKM1)=ZRNEB(JL)*ZTR1(JL)
     S              + (ZTRA1(JL,JKM1)) * ZG(JL) * (1.-ZRNEB(JL))
C
      ZREFZ(JL,1,JK)=(1.-ZRNEB(JL))*(ZRAY1(JL,JKM1)
     S                  +ZREFZ(JL,1,JKM1)*ZTRA1(JL,JKM1)*ZTRA2(JL,JKM1)
     S             /(1.-ZRAY2(JL,JKM1)*ZREFZ(JL,1,JKM1)))*ZG(JL)*ZS(JL)
     S             + ZRNEB(JL) * ZRE2(JL)
C
      ZTR(JL,1,JKM1)= ZRNEB(JL) * ZTR2(JL)
     S              + (ZTRA1(JL,JKM1)/(1.-ZRAY2(JL,JKM1)
     S              * ZREFZ(JL,1,JKM1)))
     S              * ZG(JL) * (1. -ZRNEB(JL))
C
 323  CONTINUE
 324  CONTINUE
C
C*         3.3  REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL
C               -------------------------------------------------
C
 330  CONTINUE
C
      DO 351 JREF=1,2
C
      JN = JN + 1
C
      DO 331 JL = 1, KDLON
      ZRJ(JL,JN,KFLEV+1) = 1.
      ZRK(JL,JN,KFLEV+1) = ZREFZ(JL,JREF,KFLEV+1)
 331  CONTINUE
C
      DO 333 JK = 1 , KFLEV
      JKL = KFLEV+1 - JK
      JKLP1 = JKL + 1
      DO 332 JL = 1, KDLON
      ZRE11 = ZRJ(JL,JN,JKLP1) * ZTR(JL,JREF,JKL)
      ZRJ(JL,JN,JKL) = ZRE11
      ZRK(JL,JN,JKL) = ZRE11 * ZREFZ(JL,JREF,JKL)
 332  CONTINUE
 333  CONTINUE
 351  CONTINUE
 361  CONTINUE
C
C
C     ------------------------------------------------------------------
C
C*         4.    INVERT GREY AND CONTINUUM FLUXES
C                --------------------------------
C
 400  CONTINUE
C
C
C*         4.1   UPWARD (ZRK) AND DOWNWARD (ZRJ) PSEUDO-FLUXES
C                ---------------------------------------------
C
 410  CONTINUE
C
      DO 414 JK = 1 , KFLEV+1
      DO 413 JAJ = 1 , 5 , 2
      JAJP = JAJ + 1
      DO 412 JL = 1, KDLON
      ZRJ(JL,JAJ,JK)=        ZRJ(JL,JAJ,JK) - ZRJ(JL,JAJP,JK)
      ZRK(JL,JAJ,JK)=        ZRK(JL,JAJ,JK) - ZRK(JL,JAJP,JK)
      ZRJ(JL,JAJ,JK)= MAX( ZRJ(JL,JAJ,JK) , ZEELOG )
      ZRK(JL,JAJ,JK)= MAX( ZRK(JL,JAJ,JK) , ZEELOG )
 412  CONTINUE
 413  CONTINUE
 414  CONTINUE
C
      DO 417 JK = 1 , KFLEV+1
      DO 416 JAJ = 2 , 6 , 2
      DO 415 JL = 1, KDLON
      ZRJ(JL,JAJ,JK)= MAX( ZRJ(JL,JAJ,JK) , ZEELOG )
      ZRK(JL,JAJ,JK)= MAX( ZRK(JL,JAJ,JK) , ZEELOG )
 415  CONTINUE
 416  CONTINUE
 417  CONTINUE
C
C*         4.2    EFFECTIVE ABSORBER AMOUNTS BY INVERSE LAPLACE
C                 ---------------------------------------------
C
 420  CONTINUE
C
      DO 437 JK = 1 , KFLEV+1
      JKKI = 1
      DO 425 JAJ = 1 , 2
      IIND2(1)=JAJ
      IIND2(2)=JAJ
      DO 424 JN = 1 , 2
      JN2J = JN + 2 * JAJ
      JKKP4 = JKKI + 4
C
C*         4.2.1  EFFECTIVE ABSORBER AMOUNTS
C                 --------------------------
C
 4210 CONTINUE
C
      DO 4211 JL = 1, KDLON
      ZW2(JL,1) = LOG( ZRJ(JL,JN,JK) / ZRJ(JL,JN2J,JK))
     S                               / PAKI(JL,JAJ)
      ZW2(JL,2) = LOG( ZRK(JL,JN,JK) / ZRK(JL,JN2J,JK))
     S                               / PAKI(JL,JAJ)
 4211 CONTINUE
C
C*         4.2.2  TRANSMISSION FUNCTION
C                 ---------------------
C
 4220 CONTINUE
C
      CALL SWTT1(KNU, 2, IIND2, ZW2, ZR2)
C
      DO 4221 JL = 1, KDLON
      ZRL(JL,JKKI) = ZR2(JL,1)
      ZRUEF(JL,JKKI) = ZW2(JL,1)
      ZRL(JL,JKKP4) = ZR2(JL,2)
      ZRUEF(JL,JKKP4) = ZW2(JL,2)
 4221 CONTINUE
C
      JKKI=JKKI+1
 424  CONTINUE
 425  CONTINUE
C
C*         4.3    UPWARD AND DOWNWARD FLUXES WITH H2O AND UMG ABSORPTION
C                 ------------------------------------------------------
C
 430  CONTINUE
C
      DO 431 JL = 1, KDLON
      PFDOWN(JL,JK) = ZRJ(JL,1,JK) * ZRL(JL,1) * ZRL(JL,3)
     S              + ZRJ(JL,2,JK) * ZRL(JL,2) * ZRL(JL,4)
      PFUP(JL,JK)   = ZRK(JL,1,JK) * ZRL(JL,5) * ZRL(JL,7)
     S              + ZRK(JL,2,JK) * ZRL(JL,6) * ZRL(JL,8)
 431  CONTINUE
 437  CONTINUE
C
C
C     ------------------------------------------------------------------
C
C*         5.    MOLECULAR ABSORPTION ON CLEAR-SKY FLUXES
C                ----------------------------------------
C
 500  CONTINUE
C
C
C*         5.1   DOWNWARD FLUXES
C                ---------------
C
 510  CONTINUE
C
      JAJ = 2
      IIND3(1)=1
      IIND3(2)=2
      IIND3(3)=3
C      
      DO 511 JL = 1, KDLON
      ZW3(JL,1)=0.
      ZW3(JL,2)=0.
      ZW3(JL,3)=0.
      ZW4(JL)  =0.
      ZW5(JL)  =0.
      ZR4(JL)  =1.
      ZFD(JL,KFLEV+1)= ZRJ0(JL,JAJ,KFLEV+1)
 511  CONTINUE
      DO 514 JK = 1 , KFLEV
      IKL = KFLEV+1-JK
      DO 512 JL = 1, KDLON
      ZW3(JL,1)=ZW3(JL,1)+PUD(JL,1,IKL)/ZRMU0(JL,IKL)
      ZW3(JL,2)=ZW3(JL,2)+PUD(JL,2,IKL)/ZRMU0(JL,IKL)
      ZW3(JL,3)=ZW3(JL,3)+POZ(JL,  IKL)/ZRMU0(JL,IKL)
      ZW4(JL)  =ZW4(JL)  +PUD(JL,4,IKL)/ZRMU0(JL,IKL)
      ZW5(JL)  =ZW5(JL)  +PUD(JL,5,IKL)/ZRMU0(JL,IKL)
 512  CONTINUE
C
      CALL SWTT1(KNU, 3, IIND3, ZW3, ZR3)
C
      DO 513 JL = 1, KDLON
C     ZR4(JL) = EXP(-RSWCE*ZW4(JL)-RSWCP*ZW5(JL))
      ZFD(JL,IKL) = ZR3(JL,1)*ZR3(JL,2)*ZR3(JL,3)*ZR4(JL)
     S            * ZRJ0(JL,JAJ,IKL)
 513  CONTINUE
 514  CONTINUE
C
C
C*         5.2   UPWARD FLUXES
C                -------------
C
 520  CONTINUE
C
      DO 525 JL = 1, KDLON
      ZFU(JL,1) = ZFD(JL,1)*PALBP(JL,KNU)
 525  CONTINUE
C
      DO 528 JK = 2 , KFLEV+1
      IKM1=JK-1
      DO 526 JL = 1, KDLON
      ZW3(JL,1)=ZW3(JL,1)+PUD(JL,1,IKM1)*1.66
      ZW3(JL,2)=ZW3(JL,2)+PUD(JL,2,IKM1)*1.66
      ZW3(JL,3)=ZW3(JL,3)+POZ(JL,  IKM1)*1.66
      ZW4(JL)  =ZW4(JL)  +PUD(JL,4,IKM1)*1.66
      ZW5(JL)  =ZW5(JL)  +PUD(JL,5,IKM1)*1.66
 526  CONTINUE
C
      CALL SWTT1(KNU, 3, IIND3, ZW3, ZR3)
C
      DO 527 JL = 1, KDLON
C     ZR4(JL) = EXP(-RSWCE*ZW4(JL)-RSWCP*ZW5(JL))
      ZFU(JL,JK) = ZR3(JL,1)*ZR3(JL,2)*ZR3(JL,3)*ZR4(JL)
     S           * ZRK0(JL,JAJ,JK)
 527  CONTINUE
 528  CONTINUE
C
C
C     ------------------------------------------------------------------
C
C*         6.     INTRODUCTION OF OZONE AND H2O CONTINUUM ABSORPTION
C                 --------------------------------------------------
C
 600  CONTINUE
      IABS=3
C
C*         6.1    DOWNWARD FLUXES
C                 ---------------
C
 610  CONTINUE
      DO 611 JL = 1, KDLON
      ZW1(JL)=0.
      ZW4(JL)=0.
      ZW5(JL)=0.
      ZR1(JL)=0.
      PFDOWN(JL,KFLEV+1) = ((1.-PCLEAR(JL))*PFDOWN(JL,KFLEV+1)
     S                   + PCLEAR(JL) * ZFD(JL,KFLEV+1)) * RSUN(KNU)
 611  CONTINUE
C
      DO 614 JK = 1 , KFLEV
      IKL=KFLEV+1-JK
      DO 612 JL = 1, KDLON
      ZW1(JL) = ZW1(JL)+POZ(JL,  IKL)/ZRMUE(JL,IKL)
      ZW4(JL) = ZW4(JL)+PUD(JL,4,IKL)/ZRMUE(JL,IKL)
      ZW5(JL) = ZW5(JL)+PUD(JL,5,IKL)/ZRMUE(JL,IKL)
C     ZR4(JL) = EXP(-RSWCE*ZW4(JL)-RSWCP*ZW5(JL))
 612  CONTINUE
C
      CALL SWTT(KNU, IABS, ZW1, ZR1)
C
      DO 613 JL = 1, KDLON
      PFDOWN(JL,IKL) = ((1.-PCLEAR(JL))*ZR1(JL)*ZR4(JL)*PFDOWN(JL,IKL)
     S                     +PCLEAR(JL)*ZFD(JL,IKL)) * RSUN(KNU)
 613  CONTINUE
 614  CONTINUE
C
C
C*         6.2    UPWARD FLUXES
C                 -------------
C
 620  CONTINUE
      DO 621 JL = 1, KDLON
      PFUP(JL,1) = ((1.-PCLEAR(JL))*ZR1(JL)*ZR4(JL) * PFUP(JL,1)
     S                 +PCLEAR(JL)*ZFU(JL,1)) * RSUN(KNU)
 621  CONTINUE
C
      DO 624 JK = 2 , KFLEV+1
      IKM1=JK-1
      DO 622 JL = 1, KDLON
      ZW1(JL) = ZW1(JL)+POZ(JL  ,IKM1)*1.66
      ZW4(JL) = ZW4(JL)+PUD(JL,4,IKM1)*1.66
      ZW5(JL) = ZW5(JL)+PUD(JL,5,IKM1)*1.66
C     ZR4(JL) = EXP(-RSWCE*ZW4(JL)-RSWCP*ZW5(JL))
 622  CONTINUE
C
      CALL SWTT(KNU, IABS, ZW1, ZR1)
C
      DO 623 JL = 1, KDLON
      PFUP(JL,JK) = ((1.-PCLEAR(JL))*ZR1(JL)*ZR4(JL) * PFUP(JL,JK)
     S                 +PCLEAR(JL)*ZFU(JL,JK)) * RSUN(KNU)
 623  CONTINUE
 624  CONTINUE
C
C     ------------------------------------------------------------------
C
      RETURN
      END
