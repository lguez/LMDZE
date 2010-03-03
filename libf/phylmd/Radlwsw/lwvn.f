      SUBROUTINE LWVN(KUAER,KTRAER
     R  , PABCU,PDBSL,PGA,PGB
     S  , PADJD,PADJU,PCNTRB,PDBDT)
      use dimens_m
      use dimphy
      use raddim
            use raddimlw
      IMPLICIT none
C
C-----------------------------------------------------------------------
C     PURPOSE.
C     --------
C           CARRIES OUT THE VERTICAL INTEGRATION ON NEARBY LAYERS
C           TO GIVE LONGWAVE FLUXES OR RADIANCES
C
C     METHOD.
C     -------
C
C          1. PERFORMS THE VERTICAL INTEGRATION CORRESPONDING TO THE
C     CONTRIBUTIONS OF THE ADJACENT LAYERS USING A GAUSSIAN QUADRATURE
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
C
C* ARGUMENTS:
C
      INTEGER KUAER,KTRAER
C
      REAL*8 PABCU(KDLON,NUA,3*KFLEV+1) ! ABSORBER AMOUNTS
      REAL*8 PDBSL(KDLON,Ninter,KFLEV*2) ! SUB-LAYER PLANCK FUNCTION GRADIENT
      REAL*8 PGA(KDLON,8,2,KFLEV) ! PADE APPROXIMANTS
      REAL*8 PGB(KDLON,8,2,KFLEV) ! PADE APPROXIMANTS
C
      REAL*8 PADJD(KDLON,KFLEV+1) ! CONTRIBUTION OF ADJACENT LAYERS
      REAL*8 PADJU(KDLON,KFLEV+1) ! CONTRIBUTION OF ADJACENT LAYERS
      REAL*8 PCNTRB(KDLON,KFLEV+1,KFLEV+1) ! CLEAR-SKY ENERGY EXCHANGE MATRIX
      REAL*8 PDBDT(KDLON,Ninter,KFLEV) !  LAYER PLANCK FUNCTION GRADIENT
C
C* LOCAL ARRAYS:
C
      REAL*8 ZGLAYD(KDLON)
      REAL*8 ZGLAYU(KDLON)
      REAL*8 ZTT(KDLON,NTRA)
      REAL*8 ZTT1(KDLON,NTRA)
      REAL*8 ZTT2(KDLON,NTRA)
      REAL*8 ZUU(KDLON,NUA)
C
      INTEGER jk, jl, ja, im12, ind, inu, ixu, jg
      INTEGER ixd, ibs, idd, imu, jk1, jk2, jnu
      REAL*8 zwtr
c
C* Data Block:
c
      REAL*8 WG1(2)
      SAVE WG1
      DATA (WG1(jk),jk=1,2) /1.0, 1.0/
C-----------------------------------------------------------------------
C
C*         1.    INITIALIZATION
C                --------------
C
 100  CONTINUE
C
C*         1.1     INITIALIZE LAYER CONTRIBUTIONS
C                  ------------------------------
C
 110  CONTINUE
C
      DO 112 JK = 1 , KFLEV+1
      DO 111 JL = 1, KDLON
      PADJD(JL,JK) = 0.
      PADJU(JL,JK) = 0.
 111  CONTINUE
 112  CONTINUE
C
C*         1.2     INITIALIZE TRANSMISSION FUNCTIONS
C                  ---------------------------------
C
 120  CONTINUE
C
      DO 122 JA = 1 , NTRA
      DO 121 JL = 1, KDLON
      ZTT (JL,JA) = 1.0
      ZTT1(JL,JA) = 1.0
      ZTT2(JL,JA) = 1.0
 121  CONTINUE
 122  CONTINUE
C
      DO 124 JA = 1 , NUA
      DO 123 JL = 1, KDLON
      ZUU(JL,JA) = 0.
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
C
C*         2.1     CONTRIBUTION FROM ADJACENT LAYERS
C                  ---------------------------------
C
 210  CONTINUE
C
      DO 215 JK = 1 , KFLEV
C
C*         2.1.1   DOWNWARD LAYERS
C                  ---------------
C
 2110 CONTINUE
C
      IM12 = 2 * (JK - 1)
      IND = (JK - 1) * NG1P1 + 1
      IXD = IND
      INU = JK * NG1P1 + 1
      IXU = IND
C
      DO 2111 JL = 1, KDLON
      ZGLAYD(JL) = 0.
      ZGLAYU(JL) = 0.
 2111 CONTINUE
C
      DO 213 JG = 1 , NG1
      IBS = IM12 + JG
      IDD = IXD + JG
      DO 2113 JA = 1 , KUAER
      DO 2112 JL = 1, KDLON
      ZUU(JL,JA) = PABCU(JL,JA,IND) - PABCU(JL,JA,IDD)
 2112 CONTINUE
 2113 CONTINUE
C
C
      CALL LWTT(PGA(1,1,1,JK), PGB(1,1,1,JK), ZUU, ZTT)
C
      DO 2114 JL = 1, KDLON
      ZWTR=PDBSL(JL,1,IBS)*ZTT(JL,1)          *ZTT(JL,10)
     S    +PDBSL(JL,2,IBS)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11)
     S    +PDBSL(JL,3,IBS)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12)
     S    +PDBSL(JL,4,IBS)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13)
     S    +PDBSL(JL,5,IBS)*ZTT(JL,3)          *ZTT(JL,14)
     S    +PDBSL(JL,6,IBS)*ZTT(JL,6)          *ZTT(JL,15)
      ZGLAYD(JL)=ZGLAYD(JL)+ZWTR*WG1(JG)
 2114 CONTINUE
C
C*         2.1.2   DOWNWARD LAYERS
C                  ---------------
C
 2120 CONTINUE
C
      IMU = IXU + JG
      DO 2122 JA = 1 , KUAER
      DO 2121 JL = 1, KDLON
      ZUU(JL,JA) = PABCU(JL,JA,IMU) - PABCU(JL,JA,INU)
 2121 CONTINUE
 2122 CONTINUE
C
C
      CALL LWTT(PGA(1,1,1,JK), PGB(1,1,1,JK), ZUU, ZTT)
C
      DO 2123 JL = 1, KDLON
      ZWTR=PDBSL(JL,1,IBS)*ZTT(JL,1)          *ZTT(JL,10)
     S    +PDBSL(JL,2,IBS)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11)
     S    +PDBSL(JL,3,IBS)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12)
     S    +PDBSL(JL,4,IBS)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13)
     S    +PDBSL(JL,5,IBS)*ZTT(JL,3)          *ZTT(JL,14)
     S    +PDBSL(JL,6,IBS)*ZTT(JL,6)          *ZTT(JL,15)
      ZGLAYU(JL)=ZGLAYU(JL)+ZWTR*WG1(JG)
 2123 CONTINUE
C
 213  CONTINUE
C
      DO 214 JL = 1, KDLON
      PADJD(JL,JK) = ZGLAYD(JL)
      PCNTRB(JL,JK,JK+1) = ZGLAYD(JL)
      PADJU(JL,JK+1) = ZGLAYU(JL)
      PCNTRB(JL,JK+1,JK) = ZGLAYU(JL)
      PCNTRB(JL,JK  ,JK) = 0.0
 214  CONTINUE
C
 215  CONTINUE
C
      DO 218 JK = 1 , KFLEV
      JK2 = 2 * JK
      JK1 = JK2 - 1
      DO 217 JNU = 1 , Ninter
      DO 216 JL = 1, KDLON
      PDBDT(JL,JNU,JK) = PDBSL(JL,JNU,JK1) + PDBSL(JL,JNU,JK2)
 216  CONTINUE
 217  CONTINUE
 218  CONTINUE
C
      RETURN
C
      END
