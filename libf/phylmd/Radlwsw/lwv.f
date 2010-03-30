      SUBROUTINE LWV(KUAER,KTRAER, KLIM
     R  , PABCU,PB,PBINT,PBSUIN,PBSUR,PBTOP,PDBSL,PEMIS,PPMB,PTAVE
     R  , PGA,PGB,PGASUR,PGBSUR,PGATOP,PGBTOP
     S  , PCNTRB,PCTS,PFLUC)
      use dimens_m
      use dimphy
      use YOMCST
      use raddim
            use raddimlw
      IMPLICIT none
C
C-----------------------------------------------------------------------
C     PURPOSE.
C     --------
C           CARRIES OUT THE VERTICAL INTEGRATION TO GIVE LONGWAVE
C           FLUXES OR RADIANCES
C
C     METHOD.
C     -------
C
C          1. PERFORMS THE VERTICAL INTEGRATION DISTINGUISHING BETWEEN
C     CONTRIBUTIONS BY -  THE NEARBY LAYERS
C                      -  THE DISTANT LAYERS
C                      -  THE BOUNDARY TERMS
C          2. COMPUTES THE CLEAR-SKY DOWNWARD AND UPWARD EMISSIVITIES.
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
      INTEGER KUAER,KTRAER, KLIM
C
      REAL*8 PABCU(KDLON,NUA,3*KFLEV+1) ! EFFECTIVE ABSORBER AMOUNTS
      REAL*8 PB(KDLON,Ninter,KFLEV+1) ! SPECTRAL HALF-LEVEL PLANCK FUNCTIONS
      REAL*8 PBINT(KDLON,KFLEV+1) ! HALF-LEVEL PLANCK FUNCTIONS
      REAL*8 PBSUR(KDLON,Ninter) ! SURFACE SPECTRAL PLANCK FUNCTION
      REAL*8 PBSUIN(KDLON) ! SURFACE PLANCK FUNCTION
      REAL*8 PBTOP(KDLON,Ninter) ! T.O.A. SPECTRAL PLANCK FUNCTION
      REAL*8 PDBSL(KDLON,Ninter,KFLEV*2) ! SUB-LAYER PLANCK FUNCTION GRADIENT
      REAL*8 PEMIS(KDLON) ! SURFACE EMISSIVITY
      REAL*8 PPMB(KDLON,KFLEV+1) ! HALF-LEVEL PRESSURE (MB)
      REAL*8 PTAVE(KDLON,KFLEV) ! TEMPERATURE
      REAL*8 PGA(KDLON,8,2,KFLEV) ! PADE APPROXIMANTS
      REAL*8 PGB(KDLON,8,2,KFLEV) ! PADE APPROXIMANTS
      REAL*8 PGASUR(KDLON,8,2) ! PADE APPROXIMANTS
      REAL*8 PGBSUR(KDLON,8,2) ! PADE APPROXIMANTS
      REAL*8 PGATOP(KDLON,8,2) ! PADE APPROXIMANTS
      REAL*8 PGBTOP(KDLON,8,2) ! PADE APPROXIMANTS
C
      REAL*8 PCNTRB(KDLON,KFLEV+1,KFLEV+1) ! CLEAR-SKY ENERGY EXCHANGE MATRIX
      REAL*8 PCTS(KDLON,KFLEV) ! COOLING-TO-SPACE TERM
      REAL*8 PFLUC(KDLON,2,KFLEV+1) ! CLEAR-SKY RADIATIVE FLUXES
C-----------------------------------------------------------------------
C LOCAL VARIABLES:
      REAL*8 ZADJD(KDLON,KFLEV+1)
      REAL*8 ZADJU(KDLON,KFLEV+1)
      REAL*8 ZDBDT(KDLON,Ninter,KFLEV)
      REAL*8 ZDISD(KDLON,KFLEV+1)
      REAL*8 ZDISU(KDLON,KFLEV+1)
C
      INTEGER jk, jl
C-----------------------------------------------------------------------
C
      DO 112 JK=1,KFLEV+1
      DO 111 JL=1, KDLON
      ZADJD(JL,JK)=0.
      ZADJU(JL,JK)=0.
      ZDISD(JL,JK)=0.
      ZDISU(JL,JK)=0.
 111  CONTINUE
 112  CONTINUE
C
      DO 114 JK=1,KFLEV
      DO 113 JL=1, KDLON
      PCTS(JL,JK)=0.
 113  CONTINUE
 114  CONTINUE
C
C* CONTRIBUTION FROM ADJACENT LAYERS
C
      CALL LWVN(KUAER,KTRAER
     R  , PABCU,PDBSL,PGA,PGB
     S  , ZADJD,ZADJU,PCNTRB,ZDBDT)
C* CONTRIBUTION FROM DISTANT LAYERS
C
      CALL LWVD(KUAER,KTRAER
     R  , PABCU,ZDBDT,PGA,PGB
     S  , PCNTRB,ZDISD,ZDISU)
C
C* EXCHANGE WITH THE BOUNDARIES
C
      CALL LWVB(KUAER,KTRAER, KLIM
     R  , PABCU,ZADJD,ZADJU,PB,PBINT,PBSUIN,PBSUR,PBTOP
     R  , ZDISD,ZDISU,PEMIS,PPMB
     R  , PGA,PGB,PGASUR,PGBSUR,PGATOP,PGBTOP
     S  , PCTS,PFLUC)
C
C
      RETURN
      END