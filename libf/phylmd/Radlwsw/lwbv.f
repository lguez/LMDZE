      SUBROUTINE LWBV(KLIM,PDP,PDT0,PEMIS,PPMB,PTL,PTAVE,PABCU,
     S                PFLUC,PBINT,PBSUI,PCTS,PCNTRB)
      use dimens_m
      use dimphy
      use SUPHEC_M
      use raddim
            use raddimlw
      IMPLICIT none
C
C     PURPOSE.
C     --------
C           TO COMPUTE THE PLANCK FUNCTION AND PERFORM THE
C           VERTICAL INTEGRATION. SPLIT OUT FROM LW FOR MEMORY
C           SAVING
C
C     METHOD.
C     -------
C
C          1. COMPUTES THE PLANCK FUNCTIONS ON THE INTERFACES AND THE
C     GRADIENT OF PLANCK FUNCTIONS IN THE LAYERS.
C          2. PERFORMS THE VERTICAL INTEGRATION DISTINGUISHING THE CON-
C     TRIBUTIONS OF THE ADJACENT AND DISTANT LAYERS AND THOSE FROM THE
C     BOUNDARIES.
C          3. COMPUTES THE CLEAR-SKY COOLING RATES.
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
C        MODIFICATION : 93-10-15 M.HAMRUD (SPLIT OUT FROM LW TO SAVE
C                                          MEMORY)
C-----------------------------------------------------------------------
C* ARGUMENTS:
      INTEGER KLIM
C
      REAL*8 PDP(KDLON,KFLEV)
      REAL*8 PDT0(KDLON)
      REAL*8 PEMIS(KDLON)
      REAL*8 PPMB(KDLON,KFLEV+1)
      REAL*8 PTL(KDLON,KFLEV+1)
      REAL*8 PTAVE(KDLON,KFLEV)
C
      REAL*8 PFLUC(KDLON,2,KFLEV+1)
C     
      REAL*8 PABCU(KDLON,NUA,3*KFLEV+1)
      REAL*8 PBINT(KDLON,KFLEV+1)
      REAL*8 PBSUI(KDLON)
      REAL*8 PCTS(KDLON,KFLEV)
      REAL*8 PCNTRB(KDLON,KFLEV+1,KFLEV+1)
C
C-------------------------------------------------------------------------
C
C* LOCAL VARIABLES:
      REAL*8 ZB(KDLON,Ninter,KFLEV+1)
      REAL*8 ZBSUR(KDLON,Ninter)
      REAL*8 ZBTOP(KDLON,Ninter)
      REAL*8 ZDBSL(KDLON,Ninter,KFLEV*2)
      REAL*8 ZGA(KDLON,8,2,KFLEV)
      REAL*8 ZGB(KDLON,8,2,KFLEV)
      REAL*8 ZGASUR(KDLON,8,2)
      REAL*8 ZGBSUR(KDLON,8,2)
      REAL*8 ZGATOP(KDLON,8,2)
      REAL*8 ZGBTOP(KDLON,8,2)
C
      INTEGER nuaer, ntraer
C     ------------------------------------------------------------------
C* COMPUTES PLANCK FUNCTIONS:
       CALL LWB(PDT0,PTAVE,PTL,
     S          ZB,PBINT,PBSUI,ZBSUR,ZBTOP,ZDBSL,
     S          ZGA,ZGB,ZGASUR,ZGBSUR,ZGATOP,ZGBTOP)
C     ------------------------------------------------------------------
C* PERFORMS THE VERTICAL INTEGRATION:
      NUAER = NUA
      NTRAER = NTRA
      CALL LWV(NUAER,NTRAER, KLIM
     R  , PABCU,ZB,PBINT,PBSUI,ZBSUR,ZBTOP,ZDBSL,PEMIS,PPMB,PTAVE
     R  , ZGA,ZGB,ZGASUR,ZGBSUR,ZGATOP,ZGBTOP
     S  , PCNTRB,PCTS,PFLUC)
C     ------------------------------------------------------------------
      RETURN
      END
