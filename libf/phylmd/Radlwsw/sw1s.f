      SUBROUTINE SW1S ( KNU
     S  ,  PAER  , flag_aer, tauae, pizae, cgae
     S  ,  PALBD , PALBP, PCG  , PCLD , PCLEAR, PCLDSW
     S  ,  PDSIG , POMEGA, POZ  , PRMU , PSEC , PTAU  , PUD  
     S  ,  PFD   , PFU)
      use dimens_m
      use dimphy
      use raddim
      IMPLICIT none
C
C     ------------------------------------------------------------------
C     PURPOSE.
C     --------
C
C          THIS ROUTINE COMPUTES THE SHORTWAVE RADIATION FLUXES IN TWO
C     SPECTRAL INTERVALS FOLLOWING FOUQUART AND BONNEL (1980).
C
C     METHOD.
C     -------
C
C          1. COMPUTES UPWARD AND DOWNWARD FLUXES CORRESPONDING TO
C     CONTINUUM SCATTERING
C          2. MULTIPLY BY OZONE TRANSMISSION FUNCTION
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
C
C* ARGUMENTS:
C
      INTEGER KNU
c-OB
      real*8 flag_aer
      real*8 tauae(kdlon,kflev,2)
      real*8 pizae(kdlon,kflev,2)
      real*8 cgae(kdlon,kflev,2)
      REAL*8 PAER(KDLON,KFLEV,5)
      REAL*8 PALBD(KDLON,2)
      REAL*8 PALBP(KDLON,2)
      REAL*8 PCG(KDLON,2,KFLEV)  
      REAL*8 PCLD(KDLON,KFLEV)
      REAL*8 PCLDSW(KDLON,KFLEV)
      REAL*8 PCLEAR(KDLON)
      REAL*8 PDSIG(KDLON,KFLEV)
      REAL*8 POMEGA(KDLON,2,KFLEV)
      REAL*8 POZ(KDLON,KFLEV)
      REAL*8 PRMU(KDLON)
      REAL*8 PSEC(KDLON)
      REAL*8 PTAU(KDLON,2,KFLEV)
      REAL*8 PUD(KDLON,5,KFLEV+1)
C
      REAL*8 PFD(KDLON,KFLEV+1)
      REAL*8 PFU(KDLON,KFLEV+1)
C
C* LOCAL VARIABLES:
C
      INTEGER IIND(4)
C      
      REAL*8 ZCGAZ(KDLON,KFLEV) 
      REAL*8 ZDIFF(KDLON)
      REAL*8 ZDIRF(KDLON)        
      REAL*8 ZPIZAZ(KDLON,KFLEV)
      REAL*8 ZRAYL(KDLON)
      REAL*8 ZRAY1(KDLON,KFLEV+1)
      REAL*8 ZRAY2(KDLON,KFLEV+1)
      REAL*8 ZREFZ(KDLON,2,KFLEV+1)
      REAL*8 ZRJ(KDLON,6,KFLEV+1)
      REAL*8 ZRJ0(KDLON,6,KFLEV+1)
      REAL*8 ZRK(KDLON,6,KFLEV+1)
      REAL*8 ZRK0(KDLON,6,KFLEV+1)
      REAL*8 ZRMUE(KDLON,KFLEV+1)
      REAL*8 ZRMU0(KDLON,KFLEV+1)
      REAL*8 ZR(KDLON,4)
      REAL*8 ZTAUAZ(KDLON,KFLEV)
      REAL*8 ZTRA1(KDLON,KFLEV+1)
      REAL*8 ZTRA2(KDLON,KFLEV+1)
      REAL*8 ZW(KDLON,4)
C
      INTEGER jl, jk, k, jaj, ikm1, ikl
c
c Prescribed Data:
c
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
C     ------------------------------------------------------------------
C
C*         1.     FIRST SPECTRAL INTERVAL (0.25-0.68 MICRON)
C                 ----------------------- ------------------
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
      ZRAYL(JL) =  RRAY(KNU,1) + PRMU(JL) * (RRAY(KNU,2) + PRMU(JL)
     S          * (RRAY(KNU,3) + PRMU(JL) * (RRAY(KNU,4) + PRMU(JL)
     S          * (RRAY(KNU,5) + PRMU(JL) *  RRAY(KNU,6)       ))))
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
     S  , PALBD ,PCG   ,PCLD  ,PDSIG ,POMEGA,ZRAYL
     S  , PSEC  ,PTAU
     S  , ZCGAZ ,ZPIZAZ,ZRAY1 ,ZRAY2 ,ZREFZ ,ZRJ  ,ZRK,ZRMUE
     S  , ZTAUAZ,ZTRA1 ,ZTRA2)
C
C
C     ------------------------------------------------------------------
C
C*         3.    OZONE ABSORPTION
C                ----------------
C
 300  CONTINUE
C
      IIND(1)=1
      IIND(2)=3
      IIND(3)=1
      IIND(4)=3
C      
C
C*         3.1   DOWNWARD FLUXES
C                ---------------
C
 310  CONTINUE
C
      JAJ = 2
C
      DO 311 JL = 1, KDLON
      ZW(JL,1)=0.
      ZW(JL,2)=0.
      ZW(JL,3)=0.
      ZW(JL,4)=0.
      PFD(JL,KFLEV+1)=((1.-PCLEAR(JL))*ZRJ(JL,JAJ,KFLEV+1)
     S     + PCLEAR(JL) *ZRJ0(JL,JAJ,KFLEV+1)) * RSUN(KNU)
 311  CONTINUE
      DO 314 JK = 1 , KFLEV
      IKL = KFLEV+1-JK
      DO 312 JL = 1, KDLON
      ZW(JL,1)=ZW(JL,1)+PUD(JL,1,IKL)/ZRMUE(JL,IKL)
      ZW(JL,2)=ZW(JL,2)+POZ(JL,  IKL)/ZRMUE(JL,IKL)
      ZW(JL,3)=ZW(JL,3)+PUD(JL,1,IKL)/ZRMU0(JL,IKL)
      ZW(JL,4)=ZW(JL,4)+POZ(JL,  IKL)/ZRMU0(JL,IKL)
 312  CONTINUE
C
      CALL SWTT1(KNU, 4, IIND, ZW, ZR)
C
      DO 313 JL = 1, KDLON
      ZDIFF(JL) = ZR(JL,1)*ZR(JL,2)*ZRJ(JL,JAJ,IKL)
      ZDIRF(JL) = ZR(JL,3)*ZR(JL,4)*ZRJ0(JL,JAJ,IKL)
      PFD(JL,IKL) = ((1.-PCLEAR(JL)) * ZDIFF(JL)
     S                  +PCLEAR(JL)  * ZDIRF(JL)) * RSUN(KNU)
 313  CONTINUE
 314  CONTINUE
C
C
C*         3.2   UPWARD FLUXES
C                -------------
C
 320  CONTINUE
C
      DO 325 JL = 1, KDLON
      PFU(JL,1) = ((1.-PCLEAR(JL))*ZDIFF(JL)*PALBD(JL,KNU)
     S               + PCLEAR(JL) *ZDIRF(JL)*PALBP(JL,KNU))
     S          * RSUN(KNU)
 325  CONTINUE
C
      DO 328 JK = 2 , KFLEV+1
      IKM1=JK-1
      DO 326 JL = 1, KDLON
      ZW(JL,1)=ZW(JL,1)+PUD(JL,1,IKM1)*1.66
      ZW(JL,2)=ZW(JL,2)+POZ(JL,  IKM1)*1.66
      ZW(JL,3)=ZW(JL,3)+PUD(JL,1,IKM1)*1.66
      ZW(JL,4)=ZW(JL,4)+POZ(JL,  IKM1)*1.66
 326  CONTINUE
C
      CALL SWTT1(KNU, 4, IIND, ZW, ZR)
C
      DO 327 JL = 1, KDLON
      ZDIFF(JL) = ZR(JL,1)*ZR(JL,2)*ZRK(JL,JAJ,JK)
      ZDIRF(JL) = ZR(JL,3)*ZR(JL,4)*ZRK0(JL,JAJ,JK)
      PFU(JL,JK) = ((1.-PCLEAR(JL)) * ZDIFF(JL)
     S                 +PCLEAR(JL)  * ZDIRF(JL)) * RSUN(KNU)
 327  CONTINUE
 328  CONTINUE
C
C     ------------------------------------------------------------------
C
      RETURN
      END
