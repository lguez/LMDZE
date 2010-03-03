      SUBROUTINE LWC(KLIM,PCLDLD,PCLDLU,PEMIS,PFLUC,
     R               PBINT,PBSUIN,PCTS,PCNTRB,
     S               PFLUX)
      use dimens_m
      use dimphy
      use raddim
      use radepsi
      use radopt
      IMPLICIT none
C
C     PURPOSE.
C     --------
C           INTRODUCES CLOUD EFFECTS ON LONGWAVE FLUXES OR
C           RADIANCES
C
C        EXPLICIT ARGUMENTS :
C        --------------------
C     ==== INPUTS ===
C PBINT  : (KDLON,0:KFLEV)     ; HALF LEVEL PLANCK FUNCTION
C PBSUIN : (KDLON)             ; SURFACE PLANCK FUNCTION
C PCLDLD : (KDLON,KFLEV)       ; DOWNWARD EFFECTIVE CLOUD FRACTION
C PCLDLU : (KDLON,KFLEV)       ; UPWARD EFFECTIVE CLOUD FRACTION
C PCNTRB : (KDLON,KFLEV+1,KFLEV+1); CLEAR-SKY ENERGY EXCHANGE
C PCTS   : (KDLON,KFLEV)       ; CLEAR-SKY LAYER COOLING-TO-SPACE
C PEMIS  : (KDLON)             ; SURFACE EMISSIVITY
C PFLUC
C     ==== OUTPUTS ===
C PFLUX(KDLON,2,KFLEV)         ; RADIATIVE FLUXES :
C                     1  ==>  UPWARD   FLUX TOTAL
C                     2  ==>  DOWNWARD FLUX TOTAL
C
C     METHOD.
C     -------
C
C          1. INITIALIZES ALL FLUXES TO CLEAR-SKY VALUES
C          2. EFFECT OF ONE OVERCAST UNITY EMISSIVITY CLOUD LAYER
C          3. EFFECT OF SEMI-TRANSPARENT, PARTIAL OR MULTI-LAYERED
C     CLOUDS
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
C        Voigt lines (loop 231 to 233)  - JJM & PhD - 01/96
C-----------------------------------------------------------------------
C* ARGUMENTS:
      INTEGER klim
      REAL*8 PFLUC(KDLON,2,KFLEV+1) ! CLEAR-SKY RADIATIVE FLUXES
      REAL*8 PBINT(KDLON,KFLEV+1)   ! HALF LEVEL PLANCK FUNCTION
      REAL*8 PBSUIN(KDLON)          ! SURFACE PLANCK FUNCTION
      REAL*8 PCNTRB(KDLON,KFLEV+1,KFLEV+1) !CLEAR-SKY ENERGY EXCHANGE
      REAL*8 PCTS(KDLON,KFLEV)      ! CLEAR-SKY LAYER COOLING-TO-SPACE
c
      REAL*8 PCLDLD(KDLON,KFLEV)
      REAL*8 PCLDLU(KDLON,KFLEV)
      REAL*8 PEMIS(KDLON)
C
      REAL*8 PFLUX(KDLON,2,KFLEV+1)
C-----------------------------------------------------------------------
C* LOCAL VARIABLES:
      INTEGER IMX(KDLON), IMXP(KDLON)
C
      REAL*8 ZCLEAR(KDLON),ZCLOUD(KDLON),ZDNF(KDLON,KFLEV+1,KFLEV+1)
     S  , ZFD(KDLON), ZFN10(KDLON), ZFU(KDLON)
     S  , ZUPF(KDLON,KFLEV+1,KFLEV+1)
      REAL*8 ZCLM(KDLON,KFLEV+1,KFLEV+1)
C
      INTEGER jk, jl, imaxc, imx1, imx2, jkj, jkp1, jkm1
      INTEGER jk1, jk2, jkc, jkcp1, jcloud
      INTEGER imxm1, imxp1
      REAL*8 zcfrac
C     ------------------------------------------------------------------
C
C*         1.     INITIALIZATION
C                 --------------
C
 100  CONTINUE
C
      IMAXC = 0
C
      DO 101 JL = 1, KDLON
      IMX(JL)=0
      IMXP(JL)=0
      ZCLOUD(JL) = 0.
 101  CONTINUE
C
C*         1.1    SEARCH THE LAYER INDEX OF THE HIGHEST CLOUD
C                 -------------------------------------------
C
 110  CONTINUE
C
      DO 112 JK = 1 , KFLEV
      DO 111 JL = 1, KDLON
      IMX1=IMX(JL)
      IMX2=JK
      IF (PCLDLU(JL,JK).GT.ZEPSC) THEN
         IMXP(JL)=IMX2
      ELSE
         IMXP(JL)=IMX1
      END IF
      IMAXC=MAX(IMXP(JL),IMAXC)
      IMX(JL)=IMXP(JL)
 111  CONTINUE
 112  CONTINUE
CGM*******
      IMAXC=KFLEV
CGM*******
C
      DO 114 JK = 1 , KFLEV+1
      DO 113 JL = 1, KDLON
      PFLUX(JL,1,JK) = PFLUC(JL,1,JK)
      PFLUX(JL,2,JK) = PFLUC(JL,2,JK)
 113  CONTINUE
 114  CONTINUE
C
C     ------------------------------------------------------------------
C
C*         2.      EFFECT OF CLOUDINESS ON LONGWAVE FLUXES
C                  ---------------------------------------
C
      IF (IMAXC.GT.0) THEN
C
         IMXP1 = IMAXC + 1
         IMXM1 = IMAXC - 1
C
C*         2.0     INITIALIZE TO CLEAR-SKY FLUXES
C                  ------------------------------
C
 200  CONTINUE
C
         DO 203 JK1=1,KFLEV+1
         DO 202 JK2=1,KFLEV+1
         DO 201 JL = 1, KDLON
         ZUPF(JL,JK2,JK1)=PFLUC(JL,1,JK1)
         ZDNF(JL,JK2,JK1)=PFLUC(JL,2,JK1)
 201     CONTINUE
 202     CONTINUE
 203     CONTINUE
C
C*         2.1     FLUXES FOR ONE OVERCAST UNITY EMISSIVITY CLOUD
C                  ----------------------------------------------
C
 210  CONTINUE
C
         DO 213 JKC = 1 , IMAXC
         JCLOUD=JKC
         JKCP1=JCLOUD+1
C
C*         2.1.1   ABOVE THE CLOUD
C                  ---------------
C
 2110 CONTINUE
C
         DO 2115 JK=JKCP1,KFLEV+1
         JKM1=JK-1
         DO 2111 JL = 1, KDLON
         ZFU(JL)=0.
 2111    CONTINUE
         IF (JK .GT. JKCP1) THEN
            DO 2113 JKJ=JKCP1,JKM1
            DO 2112 JL = 1, KDLON
            ZFU(JL) = ZFU(JL) + PCNTRB(JL,JK,JKJ)
 2112       CONTINUE
 2113       CONTINUE
         END IF
C
         DO 2114 JL = 1, KDLON
         ZUPF(JL,JKCP1,JK)=PBINT(JL,JK)-ZFU(JL)
 2114    CONTINUE
 2115    CONTINUE
C
C*         2.1.2   BELOW THE CLOUD
C                  ---------------
C
 2120 CONTINUE
C
         DO 2125 JK=1,JCLOUD
         JKP1=JK+1
         DO 2121 JL = 1, KDLON
         ZFD(JL)=0.
 2121    CONTINUE
C
         IF (JK .LT. JCLOUD) THEN
            DO 2123 JKJ=JKP1,JCLOUD
            DO 2122 JL = 1, KDLON
            ZFD(JL) = ZFD(JL) + PCNTRB(JL,JK,JKJ)
 2122       CONTINUE
 2123       CONTINUE
         END IF
         DO 2124 JL = 1, KDLON
         ZDNF(JL,JKCP1,JK)=-PBINT(JL,JK)-ZFD(JL)
 2124    CONTINUE
 2125    CONTINUE
C
 213     CONTINUE
C
C
C*         2.2     CLOUD COVER MATRIX
C                  ------------------
C
C*    ZCLM(JK1,JK2) IS THE OBSCURATION FACTOR BY CLOUD LAYERS BETWEEN
C     HALF-LEVELS JK1 AND JK2 AS SEEN FROM JK1
C
 220  CONTINUE
C
      DO 223 JK1 = 1 , KFLEV+1
      DO 222 JK2 = 1 , KFLEV+1
      DO 221 JL = 1, KDLON
      ZCLM(JL,JK1,JK2) = 0.
 221  CONTINUE
 222  CONTINUE
 223  CONTINUE
C
C
C
C*         2.4     CLOUD COVER BELOW THE LEVEL OF CALCULATION
C                  ------------------------------------------
C
 240  CONTINUE
C
      DO 244 JK1 = 2 , KFLEV+1
      DO 241 JL = 1, KDLON
      ZCLEAR(JL)=1.
      ZCLOUD(JL)=0.
 241  CONTINUE
      DO 243 JK = JK1 - 1 , 1 , -1
      DO 242 JL = 1, KDLON
      IF (NOVLP.EQ.1) THEN
c* maximum-random       
         ZCLEAR(JL)=ZCLEAR(JL)*(1.0-MAX(PCLDLU(JL,JK),ZCLOUD(JL)))
     *                        /(1.0-MIN(ZCLOUD(JL),1.-ZEPSEC))
         ZCLM(JL,JK1,JK) = 1.0 - ZCLEAR(JL)
         ZCLOUD(JL) = PCLDLU(JL,JK)
      ELSE IF (NOVLP.EQ.2) THEN 
c* maximum      
         ZCLOUD(JL) = MAX(ZCLOUD(JL) , PCLDLU(JL,JK))
         ZCLM(JL,JK1,JK) = ZCLOUD(JL)
      ELSE IF (NOVLP.EQ.3) THEN
c* random      
         ZCLEAR(JL) = ZCLEAR(JL)*(1.0 - PCLDLU(JL,JK))
         ZCLOUD(JL) = 1.0 - ZCLEAR(JL)
         ZCLM(JL,JK1,JK) = ZCLOUD(JL)
      END IF
 242  CONTINUE
 243  CONTINUE
 244  CONTINUE
C
C
C*         2.5     CLOUD COVER ABOVE THE LEVEL OF CALCULATION
C                  ------------------------------------------
C
 250  CONTINUE
C
      DO 254 JK1 = 1 , KFLEV
      DO 251 JL = 1, KDLON
      ZCLEAR(JL)=1.
      ZCLOUD(JL)=0.
 251  CONTINUE
      DO 253 JK = JK1 , KFLEV
      DO 252 JL = 1, KDLON
      IF (NOVLP.EQ.1) THEN
c* maximum-random       
         ZCLEAR(JL)=ZCLEAR(JL)*(1.0-MAX(PCLDLD(JL,JK),ZCLOUD(JL)))
     *                        /(1.0-MIN(ZCLOUD(JL),1.-ZEPSEC))
         ZCLM(JL,JK1,JK) = 1.0 - ZCLEAR(JL)
         ZCLOUD(JL) = PCLDLD(JL,JK)
      ELSE IF (NOVLP.EQ.2) THEN 
c* maximum      
         ZCLOUD(JL) = MAX(ZCLOUD(JL) , PCLDLD(JL,JK))
         ZCLM(JL,JK1,JK) = ZCLOUD(JL)
      ELSE IF (NOVLP.EQ.3) THEN
c* random      
         ZCLEAR(JL) = ZCLEAR(JL)*(1.0 - PCLDLD(JL,JK))
         ZCLOUD(JL) = 1.0 - ZCLEAR(JL)
         ZCLM(JL,JK1,JK) = ZCLOUD(JL)
      END IF
 252  CONTINUE
 253  CONTINUE
 254  CONTINUE
C
C
C
C*         3.      FLUXES FOR PARTIAL/MULTIPLE LAYERED CLOUDINESS
C                  ----------------------------------------------
C
 300  CONTINUE
C
C*         3.1     DOWNWARD FLUXES
C                  ---------------
C
 310  CONTINUE
C
      DO 311 JL = 1, KDLON
      PFLUX(JL,2,KFLEV+1) = 0.
 311  CONTINUE
C
      DO 317 JK1 = KFLEV , 1 , -1
C
C*                 CONTRIBUTION FROM CLEAR-SKY FRACTION
C
      DO 312 JL = 1, KDLON
      ZFD (JL) = (1. - ZCLM(JL,JK1,KFLEV)) * ZDNF(JL,1,JK1)
 312  CONTINUE
C
C*                 CONTRIBUTION FROM ADJACENT CLOUD
C
      DO 313 JL = 1, KDLON
      ZFD(JL) = ZFD(JL) + ZCLM(JL,JK1,JK1) * ZDNF(JL,JK1+1,JK1)
 313  CONTINUE
C
C*                 CONTRIBUTION FROM OTHER CLOUDY FRACTIONS
C
      DO 315 JK = KFLEV-1 , JK1 , -1
      DO 314 JL = 1, KDLON
      ZCFRAC = ZCLM(JL,JK1,JK+1) - ZCLM(JL,JK1,JK)
      ZFD(JL) =  ZFD(JL) + ZCFRAC * ZDNF(JL,JK+2,JK1)
 314  CONTINUE
 315  CONTINUE
C
      DO 316 JL = 1, KDLON
      PFLUX(JL,2,JK1) = ZFD (JL)
 316  CONTINUE
C
 317  CONTINUE
C
C
C
C
C*         3.2     UPWARD FLUX AT THE SURFACE
C                  --------------------------
C
 320  CONTINUE
C
      DO 321 JL = 1, KDLON
      PFLUX(JL,1,1) = PEMIS(JL)*PBSUIN(JL)-(1.-PEMIS(JL))*PFLUX(JL,2,1)
 321  CONTINUE
C
C
C
C*         3.3     UPWARD FLUXES
C                  -------------
C
 330  CONTINUE
C
      DO 337 JK1 = 2 , KFLEV+1
C
C*                 CONTRIBUTION FROM CLEAR-SKY FRACTION
C
      DO 332 JL = 1, KDLON
      ZFU (JL) = (1. - ZCLM(JL,JK1,1)) * ZUPF(JL,1,JK1)
 332  CONTINUE
C
C*                 CONTRIBUTION FROM ADJACENT CLOUD
C
      DO 333 JL = 1, KDLON
      ZFU(JL) =  ZFU(JL) + ZCLM(JL,JK1,JK1-1) * ZUPF(JL,JK1,JK1)
 333  CONTINUE
C
C*                 CONTRIBUTION FROM OTHER CLOUDY FRACTIONS
C
      DO 335 JK = 2 , JK1-1
      DO 334 JL = 1, KDLON
      ZCFRAC = ZCLM(JL,JK1,JK-1) - ZCLM(JL,JK1,JK)
      ZFU(JL) =  ZFU(JL) + ZCFRAC * ZUPF(JL,JK  ,JK1)
 334  CONTINUE
 335  CONTINUE
C
      DO 336 JL = 1, KDLON
      PFLUX(JL,1,JK1) = ZFU (JL)
 336  CONTINUE
C
 337  CONTINUE
C
C
      END IF
C
C
C*         2.3     END OF CLOUD EFFECT COMPUTATIONS
C
 230  CONTINUE
C
      IF (.NOT.LEVOIGT) THEN
        DO 231 JL = 1, KDLON
        ZFN10(JL) = PFLUX(JL,1,KLIM) + PFLUX(JL,2,KLIM)
 231    CONTINUE
        DO 233 JK = KLIM+1 , KFLEV+1
        DO 232 JL = 1, KDLON
        ZFN10(JL) = ZFN10(JL) + PCTS(JL,JK-1)
        PFLUX(JL,1,JK) = ZFN10(JL)
        PFLUX(JL,2,JK) = 0.0
 232    CONTINUE
 233    CONTINUE
      ENDIF
C
      RETURN
      END
