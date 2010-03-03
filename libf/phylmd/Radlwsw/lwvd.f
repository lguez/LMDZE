      SUBROUTINE LWVD(KUAER,KTRAER
     S  , PABCU,PDBDT
     R  , PGA,PGB
     S  , PCNTRB,PDISD,PDISU)
      use dimens_m
      use dimphy
      use raddim
            use raddimlw
      IMPLICIT none
C
C-----------------------------------------------------------------------
C     PURPOSE.
C     --------
C           CARRIES OUT THE VERTICAL INTEGRATION ON THE DISTANT LAYERS
C
C     METHOD.
C     -------
C
C          1. PERFORMS THE VERTICAL INTEGRATION CORRESPONDING TO THE
C     CONTRIBUTIONS OF THE DISTANT LAYERS USING TRAPEZOIDAL RULE
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
C* ARGUMENTS:
C
      INTEGER KUAER,KTRAER
C
      REAL*8 PABCU(KDLON,NUA,3*KFLEV+1) ! ABSORBER AMOUNTS
      REAL*8 PDBDT(KDLON,Ninter,KFLEV) ! LAYER PLANCK FUNCTION GRADIENT
      REAL*8 PGA(KDLON,8,2,KFLEV) ! PADE APPROXIMANTS
      REAL*8 PGB(KDLON,8,2,KFLEV) ! PADE APPROXIMANTS
C
      REAL*8 PCNTRB(KDLON,KFLEV+1,KFLEV+1) ! ENERGY EXCHANGE MATRIX
      REAL*8 PDISD(KDLON,KFLEV+1) !  CONTRIBUTION BY DISTANT LAYERS
      REAL*8 PDISU(KDLON,KFLEV+1) !  CONTRIBUTION BY DISTANT LAYERS
C
C* LOCAL VARIABLES:
C
      REAL*8 ZGLAYD(KDLON)
      REAL*8 ZGLAYU(KDLON)
      REAL*8 ZTT(KDLON,NTRA)
      REAL*8 ZTT1(KDLON,NTRA)
      REAL*8 ZTT2(KDLON,NTRA)
C
      INTEGER jl, jk, ja, ikp1, ikn, ikd1, jkj, ikd2
      INTEGER ikjp1, ikm1, ikj, jlk, iku1, ijkl, iku2
      INTEGER ind1, ind2, ind3, ind4, itt
      REAL*8 zww, zdzxdg, zdzxmg
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
      DO 112 JK = 1, KFLEV+1
      DO 111 JL = 1, KDLON
      PDISD(JL,JK) = 0.
      PDISU(JL,JK) = 0.
  111 CONTINUE
  112 CONTINUE
C
C*         1.2     INITIALIZE TRANSMISSION FUNCTIONS
C                  ---------------------------------
C
 120  CONTINUE
C
C
      DO 122 JA = 1, NTRA
      DO 121 JL = 1, KDLON
      ZTT (JL,JA) = 1.0
      ZTT1(JL,JA) = 1.0
      ZTT2(JL,JA) = 1.0
  121 CONTINUE
  122 CONTINUE
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
C*         2.2     CONTRIBUTION FROM DISTANT LAYERS
C                  ---------------------------------
C
 220  CONTINUE
C
C
C*         2.2.1   DISTANT AND ABOVE LAYERS
C                  ------------------------
C
 2210 CONTINUE
C
C
C
C*         2.2.2   FIRST UPPER LEVEL
C                  -----------------
C
 2220 CONTINUE
C
      DO 225 JK = 1 , KFLEV-1
      IKP1=JK+1
      IKN=(JK-1)*NG1P1+1
      IKD1= JK  *NG1P1+1
C
      CALL LWTTM(PGA(1,1,1,JK), PGB(1,1,1,JK)
     2          , PABCU(1,1,IKN),PABCU(1,1,IKD1),ZTT1)
C
C
C
C*         2.2.3   HIGHER UP
C                  ---------
C
 2230 CONTINUE
C
      ITT=1
      DO 224 JKJ=IKP1,KFLEV
      IF(ITT.EQ.1) THEN
         ITT=2
      ELSE
         ITT=1
      ENDIF
      IKJP1=JKJ+1
      IKD2= JKJ  *NG1P1+1
C
      IF(ITT.EQ.1) THEN
         CALL LWTTM(PGA(1,1,1,JKJ),PGB(1,1,1,JKJ)
     2             , PABCU(1,1,IKN),PABCU(1,1,IKD2),ZTT1)
      ELSE
         CALL LWTTM(PGA(1,1,1,JKJ),PGB(1,1,1,JKJ)
     2             , PABCU(1,1,IKN),PABCU(1,1,IKD2),ZTT2)
      ENDIF
C
      DO 2235 JA = 1, KTRAER
      DO 2234 JL = 1, KDLON
      ZTT(JL,JA) = (ZTT1(JL,JA)+ZTT2(JL,JA))*0.5
 2234 CONTINUE
 2235 CONTINUE
C
      DO 2236 JL = 1, KDLON
      ZWW=PDBDT(JL,1,JKJ)*ZTT(JL,1)          *ZTT(JL,10)
     S   +PDBDT(JL,2,JKJ)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11)
     S   +PDBDT(JL,3,JKJ)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12)
     S   +PDBDT(JL,4,JKJ)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13)
     S   +PDBDT(JL,5,JKJ)*ZTT(JL,3)          *ZTT(JL,14)
     S   +PDBDT(JL,6,JKJ)*ZTT(JL,6)          *ZTT(JL,15)
      ZGLAYD(JL)=ZWW
      ZDZXDG=ZGLAYD(JL)
      PDISD(JL,JK)=PDISD(JL,JK)+ZDZXDG
      PCNTRB(JL,JK,IKJP1)=ZDZXDG
 2236 CONTINUE
C
C
 224  CONTINUE
 225  CONTINUE
C
C
C*         2.2.4   DISTANT AND BELOW LAYERS
C                  ------------------------
C
 2240 CONTINUE
C
C
C
C*         2.2.5   FIRST LOWER LEVEL
C                  -----------------
C
 2250 CONTINUE
C
      DO 228 JK=3,KFLEV+1
      IKN=(JK-1)*NG1P1+1
      IKM1=JK-1
      IKJ=JK-2
      IKU1= IKJ  *NG1P1+1
C
C
      CALL LWTTM(PGA(1,1,1,IKJ),PGB(1,1,1,IKJ)
     2          , PABCU(1,1,IKU1),PABCU(1,1,IKN),ZTT1)
C
C
C
C*         2.2.6   DOWN BELOW
C                  ----------
C
 2260 CONTINUE
C
      ITT=1
      DO 227 JLK=1,IKJ
      IF(ITT.EQ.1) THEN
         ITT=2
      ELSE
         ITT=1
      ENDIF
      IJKL=IKM1-JLK
      IKU2=(IJKL-1)*NG1P1+1
C
C
      IF(ITT.EQ.1) THEN
         CALL LWTTM(PGA(1,1,1,IJKL),PGB(1,1,1,IJKL)
     2             , PABCU(1,1,IKU2),PABCU(1,1,IKN),ZTT1)
      ELSE
         CALL LWTTM(PGA(1,1,1,IJKL),PGB(1,1,1,IJKL)
     2             , PABCU(1,1,IKU2),PABCU(1,1,IKN),ZTT2)
      ENDIF
C
      DO 2265 JA = 1, KTRAER
      DO 2264 JL = 1, KDLON
      ZTT(JL,JA) = (ZTT1(JL,JA)+ZTT2(JL,JA))*0.5
 2264 CONTINUE
 2265 CONTINUE
C
      DO 2266 JL = 1, KDLON
      ZWW=PDBDT(JL,1,IJKL)*ZTT(JL,1)          *ZTT(JL,10)
     S   +PDBDT(JL,2,IJKL)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11)
     S   +PDBDT(JL,3,IJKL)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12)
     S   +PDBDT(JL,4,IJKL)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13)
     S   +PDBDT(JL,5,IJKL)*ZTT(JL,3)          *ZTT(JL,14)
     S   +PDBDT(JL,6,IJKL)*ZTT(JL,6)          *ZTT(JL,15)
      ZGLAYU(JL)=ZWW
      ZDZXMG=ZGLAYU(JL)
      PDISU(JL,JK)=PDISU(JL,JK)+ZDZXMG
      PCNTRB(JL,JK,IJKL)=ZDZXMG
 2266 CONTINUE
C
C
 227  CONTINUE
 228  CONTINUE
C
      RETURN
      END
