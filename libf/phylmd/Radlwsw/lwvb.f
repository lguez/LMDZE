      SUBROUTINE LWVB(KUAER,KTRAER, KLIM
     R  , PABCU,PADJD,PADJU,PB,PBINT,PBSUI,PBSUR,PBTOP
     R  , PDISD,PDISU,PEMIS,PPMB
     R  , PGA,PGB,PGASUR,PGBSUR,PGATOP,PGBTOP
     S  , PCTS,PFLUC)
      use dimens_m
      use dimphy
      use raddim
      use radopt
            use raddimlw
      IMPLICIT none
C
C-----------------------------------------------------------------------
C     PURPOSE.
C     --------
C           INTRODUCES THE EFFECTS OF THE BOUNDARIES IN THE VERTICAL
C           INTEGRATION
C
C     METHOD.
C     -------
C
C          1. COMPUTES THE ENERGY EXCHANGE WITH TOP AND SURFACE OF THE
C     ATMOSPHERE
C          2. COMPUTES THE COOLING-TO-SPACE AND HEATING-FROM-GROUND
C     TERMS FOR THE APPROXIMATE COOLING RATE ABOVE 10 HPA
C          3. ADDS UP ALL CONTRIBUTIONS TO GET THE CLEAR-SKY FLUXES
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
C        Voigt lines (loop 2413 to 2427)  - JJM & PhD - 01/96
C-----------------------------------------------------------------------
C
C*       0.1   ARGUMENTS
C              ---------
C
      INTEGER KUAER,KTRAER, KLIM
C
      REAL*8 PABCU(KDLON,NUA,3*KFLEV+1) ! ABSORBER AMOUNTS
      REAL*8 PADJD(KDLON,KFLEV+1) ! CONTRIBUTION BY ADJACENT LAYERS
      REAL*8 PADJU(KDLON,KFLEV+1) ! CONTRIBUTION BY ADJACENT LAYERS
      REAL*8 PB(KDLON,Ninter,KFLEV+1) ! SPECTRAL HALF-LEVEL PLANCK FUNCTIONS
      REAL*8 PBINT(KDLON,KFLEV+1) ! HALF-LEVEL PLANCK FUNCTIONS
      REAL*8 PBSUR(KDLON,Ninter) ! SPECTRAL SURFACE PLANCK FUNCTION
      REAL*8 PBSUI(KDLON) ! SURFACE PLANCK FUNCTION
      REAL*8 PBTOP(KDLON,Ninter) ! SPECTRAL T.O.A. PLANCK FUNCTION
      REAL*8 PDISD(KDLON,KFLEV+1) ! CONTRIBUTION BY DISTANT LAYERS
      REAL*8 PDISU(KDLON,KFLEV+1) ! CONTRIBUTION BY DISTANT LAYERS
      REAL*8 PEMIS(KDLON) ! SURFACE EMISSIVITY
      REAL*8 PPMB(KDLON,KFLEV+1) ! PRESSURE MB
      REAL*8 PGA(KDLON,8,2,KFLEV) ! PADE APPROXIMANTS
      REAL*8 PGB(KDLON,8,2,KFLEV) ! PADE APPROXIMANTS
      REAL*8 PGASUR(KDLON,8,2) ! SURFACE PADE APPROXIMANTS
      REAL*8 PGBSUR(KDLON,8,2) ! SURFACE PADE APPROXIMANTS
      REAL*8 PGATOP(KDLON,8,2) ! T.O.A. PADE APPROXIMANTS
      REAL*8 PGBTOP(KDLON,8,2) ! T.O.A. PADE APPROXIMANTS
C
      REAL*8 PFLUC(KDLON,2,KFLEV+1) ! CLEAR-SKY RADIATIVE FLUXES
      REAL*8 PCTS(KDLON,KFLEV) ! COOLING-TO-SPACE TERM
C
C* LOCAL VARIABLES:
C
      REAL*8 ZBGND(KDLON)
      REAL*8 ZFD(KDLON)
      REAL*8  ZFN10(KDLON)
      REAL*8 ZFU(KDLON)
      REAL*8  ZTT(KDLON,NTRA)
      REAL*8 ZTT1(KDLON,NTRA)
      REAL*8 ZTT2(KDLON,NTRA)
      REAL*8  ZUU(KDLON,NUA) 
      REAL*8 ZCNSOL(KDLON)
      REAL*8 ZCNTOP(KDLON)
C
      INTEGER jk, jl, ja
      INTEGER jstra, jstru
      INTEGER ind1, ind2, ind3, ind4, in, jlim
      REAL*8 zctstr
C-----------------------------------------------------------------------
C
C*         1.    INITIALIZATION
C                --------------
C
 100  CONTINUE
C
C
C*         1.2     INITIALIZE TRANSMISSION FUNCTIONS
C                  ---------------------------------
C
 120  CONTINUE
C
      DO 122 JA=1,NTRA
      DO 121 JL=1, KDLON
      ZTT (JL,JA)=1.0
      ZTT1(JL,JA)=1.0
      ZTT2(JL,JA)=1.0
 121  CONTINUE
 122  CONTINUE
C
      DO 124 JA=1,NUA
      DO 123 JL=1, KDLON
      ZUU(JL,JA)=1.0
 123  CONTINUE
 124  CONTINUE
C
C     ------------------------------------------------------------------
C
C*         2.      VERTICAL INTEGRATION
C                  --------------------
C
 200  CONTINUE
C
      IND1=0
      IND3=0
      IND4=1
      IND2=1
C
C
C*         2.3     EXCHANGE WITH TOP OF THE ATMOSPHERE
C                  -----------------------------------
C
 230  CONTINUE
C
      DO 235 JK = 1 , KFLEV
      IN=(JK-1)*NG1P1+1
C
      DO 232 JA=1,KUAER
      DO 231 JL=1, KDLON
      ZUU(JL,JA)=PABCU(JL,JA,IN)
 231  CONTINUE
 232  CONTINUE
C
C
      CALL LWTT(PGATOP(1,1,1), PGBTOP(1,1,1), ZUU, ZTT)
C
      DO 234 JL = 1, KDLON
      ZCNTOP(JL)=PBTOP(JL,1)*ZTT(JL,1)          *ZTT(JL,10)
     2      +PBTOP(JL,2)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11)
     3      +PBTOP(JL,3)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12)
     4      +PBTOP(JL,4)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13)
     5      +PBTOP(JL,5)*ZTT(JL,3)          *ZTT(JL,14)
     6      +PBTOP(JL,6)*ZTT(JL,6)          *ZTT(JL,15)
      ZFD(JL)=ZCNTOP(JL)-PBINT(JL,JK)-PDISD(JL,JK)-PADJD(JL,JK)
      PFLUC(JL,2,JK)=ZFD(JL)
 234  CONTINUE
C
 235  CONTINUE
C
      JK = KFLEV+1
      IN=(JK-1)*NG1P1+1
C
      DO 236 JL = 1, KDLON
      ZCNTOP(JL)= PBTOP(JL,1)
     1   + PBTOP(JL,2)
     2   + PBTOP(JL,3)
     3   + PBTOP(JL,4)
     4   + PBTOP(JL,5)
     5   + PBTOP(JL,6)
      ZFD(JL)=ZCNTOP(JL)-PBINT(JL,JK)-PDISD(JL,JK)-PADJD(JL,JK)
      PFLUC(JL,2,JK)=ZFD(JL)
 236  CONTINUE
C
C*         2.4     COOLING-TO-SPACE OF LAYERS ABOVE 10 HPA
C                  ---------------------------------------
C
 240  CONTINUE
C
C
C*         2.4.1   INITIALIZATION
C                  --------------
C
 2410 CONTINUE
C
      JLIM = KFLEV
C
      IF (.NOT.LEVOIGT) THEN
      DO 2412 JK = KFLEV,1,-1
      IF(PPMB(1,JK).LT.10.0) THEN
         JLIM=JK
      ENDIF   
 2412 CONTINUE
      ENDIF
      KLIM=JLIM
C
      IF (.NOT.LEVOIGT) THEN
        DO 2414 JA=1,KTRAER
        DO 2413 JL=1, KDLON
        ZTT1(JL,JA)=1.0
 2413   CONTINUE
 2414   CONTINUE
C
C*         2.4.2   LOOP OVER LAYERS ABOVE 10 HPA
C                  -----------------------------
C
 2420   CONTINUE
C
        DO 2427 JSTRA = KFLEV,JLIM,-1
        JSTRU=(JSTRA-1)*NG1P1+1
C
        DO 2423 JA=1,KUAER
        DO 2422 JL=1, KDLON
        ZUU(JL,JA)=PABCU(JL,JA,JSTRU)
 2422   CONTINUE
 2423   CONTINUE
C
C
        CALL LWTT(PGA(1,1,1,JSTRA), PGB(1,1,1,JSTRA), ZUU, ZTT)
C
        DO 2424 JL = 1, KDLON
        ZCTSTR =
     1   (PB(JL,1,JSTRA)+PB(JL,1,JSTRA+1))
     1       *(ZTT1(JL,1)           *ZTT1(JL,10)
     1       - ZTT (JL,1)           *ZTT (JL,10))
     2  +(PB(JL,2,JSTRA)+PB(JL,2,JSTRA+1))
     2       *(ZTT1(JL,2)*ZTT1(JL,7)*ZTT1(JL,11)
     2       - ZTT (JL,2)*ZTT (JL,7)*ZTT (JL,11))
     3  +(PB(JL,3,JSTRA)+PB(JL,3,JSTRA+1))
     3       *(ZTT1(JL,4)*ZTT1(JL,8)*ZTT1(JL,12)
     3       - ZTT (JL,4)*ZTT (JL,8)*ZTT (JL,12))
     4  +(PB(JL,4,JSTRA)+PB(JL,4,JSTRA+1))
     4       *(ZTT1(JL,5)*ZTT1(JL,9)*ZTT1(JL,13)
     4       - ZTT (JL,5)*ZTT (JL,9)*ZTT (JL,13))
     5  +(PB(JL,5,JSTRA)+PB(JL,5,JSTRA+1))
     5       *(ZTT1(JL,3)           *ZTT1(JL,14)
     5       - ZTT (JL,3)           *ZTT (JL,14))
     6  +(PB(JL,6,JSTRA)+PB(JL,6,JSTRA+1))
     6       *(ZTT1(JL,6)           *ZTT1(JL,15)
     6       - ZTT (JL,6)           *ZTT (JL,15))
        PCTS(JL,JSTRA)=ZCTSTR*0.5
 2424   CONTINUE
        DO 2426 JA=1,KTRAER
        DO 2425 JL=1, KDLON
        ZTT1(JL,JA)=ZTT(JL,JA)
 2425   CONTINUE
 2426   CONTINUE
 2427   CONTINUE
      ENDIF
C Mise a zero de securite pour PCTS en cas de LEVOIGT
      IF(LEVOIGT)THEN
        DO 2429 JSTRA = 1,KFLEV
        DO 2428 JL = 1, KDLON
          PCTS(JL,JSTRA)=0.
 2428   CONTINUE
 2429   CONTINUE
      ENDIF
C
C
C*         2.5     EXCHANGE WITH LOWER LIMIT
C                  -------------------------
C
 250  CONTINUE
C
      DO 251 JL = 1, KDLON
      ZBGND(JL)=PBSUI(JL)*PEMIS(JL)-(1.-PEMIS(JL))
     S               *PFLUC(JL,2,1)-PBINT(JL,1)
 251  CONTINUE
C
      JK = 1
      IN=(JK-1)*NG1P1+1
C
      DO 252 JL = 1, KDLON
      ZCNSOL(JL)=PBSUR(JL,1)
     1 +PBSUR(JL,2)
     2 +PBSUR(JL,3)
     3 +PBSUR(JL,4)
     4 +PBSUR(JL,5)
     5 +PBSUR(JL,6)
      ZCNSOL(JL)=ZCNSOL(JL)*ZBGND(JL)/PBSUI(JL)
      ZFU(JL)=ZCNSOL(JL)+PBINT(JL,JK)-PDISU(JL,JK)-PADJU(JL,JK)
      PFLUC(JL,1,JK)=ZFU(JL)
 252  CONTINUE
C
      DO 257 JK = 2 , KFLEV+1
      IN=(JK-1)*NG1P1+1
C
C
      DO 255 JA=1,KUAER
      DO 254 JL=1, KDLON
      ZUU(JL,JA)=PABCU(JL,JA,1)-PABCU(JL,JA,IN)
 254  CONTINUE
 255  CONTINUE
C
C
      CALL LWTT(PGASUR(1,1,1), PGBSUR(1,1,1), ZUU, ZTT)
C
      DO 256 JL = 1, KDLON
      ZCNSOL(JL)=PBSUR(JL,1)*ZTT(JL,1)          *ZTT(JL,10)
     2      +PBSUR(JL,2)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11)
     3      +PBSUR(JL,3)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12)
     4      +PBSUR(JL,4)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13)
     5      +PBSUR(JL,5)*ZTT(JL,3)          *ZTT(JL,14)
     6      +PBSUR(JL,6)*ZTT(JL,6)          *ZTT(JL,15)
      ZCNSOL(JL)=ZCNSOL(JL)*ZBGND(JL)/PBSUI(JL)
      ZFU(JL)=ZCNSOL(JL)+PBINT(JL,JK)-PDISU(JL,JK)-PADJU(JL,JK)
      PFLUC(JL,1,JK)=ZFU(JL)
 256  CONTINUE
C
C
 257  CONTINUE
C
C
C
C*         2.7     CLEAR-SKY FLUXES
C                  ----------------
C
 270  CONTINUE
C
      IF (.NOT.LEVOIGT) THEN
      DO 271 JL = 1, KDLON
      ZFN10(JL) = PFLUC(JL,1,JLIM) + PFLUC(JL,2,JLIM)
 271  CONTINUE
      DO 273 JK = JLIM+1,KFLEV+1
      DO 272 JL = 1, KDLON
      ZFN10(JL) = ZFN10(JL) + PCTS(JL,JK-1)
      PFLUC(JL,1,JK) = ZFN10(JL)
      PFLUC(JL,2,JK) = 0.
 272  CONTINUE
 273  CONTINUE
      ENDIF
C
C     ------------------------------------------------------------------
C
      RETURN
      END
