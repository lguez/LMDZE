!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/advxp.F,v 1.1.1.1 2004/05/19 12:53:06 lmdzadmin Exp $
!
       SUBROUTINE ADVXP(LIMIT,DTX,PBARU,SM,S0,SSX,SY,SZ
     .                ,SSXX,SSXY,SSXZ,SYY,SYZ,SZZ,ntra)
       use dimens_m
      use paramet_m
      use comconst
      use comvert
       IMPLICIT NONE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                 C
C  second-order moments (SOM) advection of tracer in X direction  C
C                                                                 C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  parametres principaux du modele
C

       INTEGER ntra
c      PARAMETER (ntra = 1)
C
C  definition de la grille du modele
C
      REAL dtx
      REAL pbaru ( iip1,jjp1,llm )
C
C  moments: SM  total mass in each grid box
C           S0  mass of tracer in each grid box
C           Si  1rst order moment in i direction
C           Sij 2nd  order moment in i and j directions
C
      REAL SM(iip1,jjp1,llm)
     +    ,S0(iip1,jjp1,llm,ntra)
      REAL SSX(iip1,jjp1,llm,ntra)
     +    ,SY(iip1,jjp1,llm,ntra)
     +    ,SZ(iip1,jjp1,llm,ntra)
      REAL SSXX(iip1,jjp1,llm,ntra)
     +    ,SSXY(iip1,jjp1,llm,ntra)
     +    ,SSXZ(iip1,jjp1,llm,ntra)
     +    ,SYY(iip1,jjp1,llm,ntra)
     +    ,SYZ(iip1,jjp1,llm,ntra)
     +    ,SZZ(iip1,jjp1,llm,ntra)

C  Local :
C  -------

C  mass fluxes across the boundaries (UGRI,VGRI,WGRI)
C  mass fluxes in kg
C  declaration :

       REAL UGRI(iip1,jjp1,llm)

C  Rem : VGRI et WGRI ne sont pas utilises dans
C  cette subroutine ( advection en x uniquement )
C
C
C  Tij are the moments for the current latitude and level
C
      REAL TM (iim)
      REAL T0 (iim,NTRA),TX (iim,NTRA)
      REAL TY (iim,NTRA),TZ (iim,NTRA)
      REAL TXX(iim,NTRA),TXY(iim,NTRA)
      REAL TXZ(iim,NTRA),TYY(iim,NTRA)
      REAL TYZ(iim,NTRA),TZZ(iim,NTRA)
C
C  the moments F are similarly defined and used as temporary
C  storage for portions of the grid boxes in transit
C
      REAL FM (iim)
      REAL F0 (iim,NTRA),FX (iim,NTRA)
      REAL FY (iim,NTRA),FZ (iim,NTRA)
      REAL FXX(iim,NTRA),FXY(iim,NTRA)
      REAL FXZ(iim,NTRA),FYY(iim,NTRA)
      REAL FYZ(iim,NTRA),FZZ(iim,NTRA)
C
C  work arrays
C
      REAL ALF (iim),ALF1(iim),ALFQ(iim),ALF1Q(iim)
      REAL ALF2(iim),ALF3(iim),ALF4(iim)
C
      REAL SMNEW(iim),UEXT(iim)
      REAL sqi,sqf
      REAL TEMPTM
      REAL SLPMAX
      REAL S1MAX,S1NEW,S2NEW

      LOGICAL LIMIT
      INTEGER NUM(jjp1),LONK,NUMK
      INTEGER lon,lati,latf,niv
      INTEGER i,i2,i3,j,jv,l,k,iter

      lon = iim
      lati=2
      latf = jjm
      niv = llm

C *** Test de passage d'arguments ******

c      DO 399 l = 1, llm
c       DO 399 j = 1, jjp1
c        DO 399 i = 1, iip1
c         IF (S0(i,j,l,ntra) .lt. 0. ) THEN
c         PRINT*,'S0(',i,j,l,')=',S0(i,j,l,ntra)
c	     print*, 'SSX(',i,j,l,')=',SSX(i,j,l,ntra)
c         print*, 'SY(',i,j,l,')=',SY(i,j,l,ntra)
c         print*, 'SZ(',i,j,l,')=',SZ(i,j,l,ntra)
c         PRINT*, 'AIE !! debut ADVXP - pbl arg. passage dans ADVXP'
cc            STOP
c         ENDIF
c  399 CONTINUE

C *** Test : diagnostique de la qtite totale de traceur
C            dans l'atmosphere avant l'advection
c
      sqi =0.
      sqf =0.
c
      DO l = 1, llm
      DO j = 1, jjp1
      DO i = 1, iim
	 sqi = sqi + S0(i,j,l,ntra)
      END DO
      END DO
      END DO
      PRINT*,'------ DIAG DANS ADVX2 - ENTREE -----'
      PRINT*,'sqi=',sqi
c test
c  -------------------------------------
        DO 300 j =1,jjp1
         NUM(j) =1 
 300  CONTINUE
c       DO l=1,llm
c      NUM(2,l)=6
c      NUM(3,l)=6
c      NUM(jjm-1,l)=6  
c      NUM(jjm,l)=6
c      ENDDO
c        DO j=2,6
c       NUM(j)=12
c       ENDDO
c       DO j=jjm-5,jjm-1 
c       NUM(j)=12
c       ENDDO

C  Interface : adaptation nouveau modele
C  -------------------------------------
C
C  ---------------------------------------------------------
C  Conversion des flux de masses en kg/s
C  pbaru est en N/s d'ou :
C  ugri est en kg/s

       DO 500 l = 1,llm
       DO 500 j = 1,jjp1
       DO 500 i = 1,iip1
       ugri (i,j,llm+1-l) =pbaru (i,j,l) 
 500   CONTINUE

C  ---------------------------------------------------------
C  start here
C
C  boucle principale sur les niveaux et les latitudes
C     
      DO 1 L=1,NIV
      DO 1 K=lati,latf

C
C  initialisation
C
C  program assumes periodic boundaries in X
C
      DO 10 I=2,LON
         SMNEW(I)=SM(I,K,L)+(UGRI(I-1,K,L)-UGRI(I,K,L))*DTX
 10   CONTINUE
      SMNEW(1)=SM(1,K,L)+(UGRI(LON,K,L)-UGRI(1,K,L))*DTX
C
C  modifications for extended polar zones
C
      NUMK=NUM(K)
      LONK=LON/NUMK
C
      IF(NUMK.GT.1) THEN
C
      DO 111 I=1,LON
         TM(I)=0.
 111  CONTINUE
      DO 112 JV=1,NTRA
      DO 1120 I=1,LON
         T0 (I,JV)=0.
         TX (I,JV)=0.
         TY (I,JV)=0.
         TZ (I,JV)=0.
         TXX(I,JV)=0.
         TXY(I,JV)=0.
         TXZ(I,JV)=0.
         TYY(I,JV)=0.
         TYZ(I,JV)=0.
         TZZ(I,JV)=0.
 1120 CONTINUE
 112  CONTINUE
C
      DO 11 I2=1,NUMK
C
         DO 113 I=1,LONK
            I3=(I-1)*NUMK+I2
            TM(I)=TM(I)+SM(I3,K,L)
            ALF(I)=SM(I3,K,L)/TM(I)
            ALF1(I)=1.-ALF(I)
            ALFQ(I)=ALF(I)*ALF(I)
            ALF1Q(I)=ALF1(I)*ALF1(I)
            ALF2(I)=ALF1(I)-ALF(I)
            ALF3(I)=ALF(I)*ALF1(I)
 113     CONTINUE
C
         DO 114 JV=1,NTRA
         DO 1140 I=1,LONK
            I3=(I-1)*NUMK+I2
            TEMPTM=-ALF(I)*T0(I,JV)+ALF1(I)*S0(I3,K,L,JV)
            T0 (I,JV)=T0(I,JV)+S0(I3,K,L,JV)
            TXX(I,JV)=ALFQ(I)*SSXX(I3,K,L,JV)+ALF1Q(I)*TXX(I,JV)
     +        +5.*( ALF3(I)*(SSX(I3,K,L,JV)-TX(I,JV))+ALF2(I)*TEMPTM )
            TX (I,JV)=ALF(I)*SSX(I3,K,L,JV)+ALF1(I)*TX(I,JV)+3.*TEMPTM
            TXY(I,JV)=ALF (I)*SSXY(I3,K,L,JV)+ALF1(I)*TXY(I,JV)
     +           +3.*(ALF1(I)*SY (I3,K,L,JV)-ALF (I)*TY (I,JV))
            TXZ(I,JV)=ALF (I)*SSXZ(I3,K,L,JV)+ALF1(I)*TXZ(I,JV)
     +           +3.*(ALF1(I)*SZ (I3,K,L,JV)-ALF (I)*TZ (I,JV))
            TY (I,JV)=TY (I,JV)+SY (I3,K,L,JV)
            TZ (I,JV)=TZ (I,JV)+SZ (I3,K,L,JV)
            TYY(I,JV)=TYY(I,JV)+SYY(I3,K,L,JV)
            TYZ(I,JV)=TYZ(I,JV)+SYZ(I3,K,L,JV)
            TZZ(I,JV)=TZZ(I,JV)+SZZ(I3,K,L,JV)
 1140    CONTINUE
 114     CONTINUE
C
 11   CONTINUE
C
      ELSE
C
      DO 115 I=1,LON
         TM(I)=SM(I,K,L)
 115  CONTINUE
      DO 116 JV=1,NTRA
      DO 1160 I=1,LON
         T0 (I,JV)=S0 (I,K,L,JV)
         TX (I,JV)=SSX (I,K,L,JV)
         TY (I,JV)=SY (I,K,L,JV)
         TZ (I,JV)=SZ (I,K,L,JV)
         TXX(I,JV)=SSXX(I,K,L,JV)
         TXY(I,JV)=SSXY(I,K,L,JV)
         TXZ(I,JV)=SSXZ(I,K,L,JV)
         TYY(I,JV)=SYY(I,K,L,JV)
         TYZ(I,JV)=SYZ(I,K,L,JV)
         TZZ(I,JV)=SZZ(I,K,L,JV)
 1160 CONTINUE
 116  CONTINUE
C
      ENDIF
C
      DO 117 I=1,LONK
         UEXT(I)=UGRI(I*NUMK,K,L)
 117  CONTINUE
C
C  place limits on appropriate moments before transport
C      (if flux-limiting is to be applied)
C
      IF(.NOT.LIMIT) GO TO 13
C
      DO 12 JV=1,NTRA
      DO 120 I=1,LONK
        IF(T0(I,JV).GT.0.) THEN
          SLPMAX=T0(I,JV)
          S1MAX=1.5*SLPMAX
          S1NEW=AMIN1(S1MAX,AMAX1(-S1MAX,TX(I,JV)))
          S2NEW=AMIN1( 2.*SLPMAX-ABS(S1NEW)/3. ,
     +                 AMAX1(ABS(S1NEW)-SLPMAX,TXX(I,JV)) )
          TX (I,JV)=S1NEW
          TXX(I,JV)=S2NEW
          TXY(I,JV)=AMIN1(SLPMAX,AMAX1(-SLPMAX,TXY(I,JV)))
          TXZ(I,JV)=AMIN1(SLPMAX,AMAX1(-SLPMAX,TXZ(I,JV)))
        ELSE
          TX (I,JV)=0.
          TXX(I,JV)=0.
          TXY(I,JV)=0.
          TXZ(I,JV)=0.
        ENDIF
 120  CONTINUE
 12   CONTINUE
C
 13   CONTINUE
C
C  calculate flux and moments between adjacent boxes
C  1- create temporary moments/masses for partial boxes in transit
C  2- reajusts moments remaining in the box
C
C  flux from IP to I if U(I).lt.0
C
      DO 140 I=1,LONK-1
         IF(UEXT(I).LT.0.) THEN
           FM(I)=-UEXT(I)*DTX
           ALF(I)=FM(I)/TM(I+1)
           TM(I+1)=TM(I+1)-FM(I)
         ENDIF
 140  CONTINUE
C
      I=LONK
      IF(UEXT(I).LT.0.) THEN
        FM(I)=-UEXT(I)*DTX
        ALF(I)=FM(I)/TM(1)
        TM(1)=TM(1)-FM(I)
      ENDIF
C
C  flux from I to IP if U(I).gt.0
C
      DO 141 I=1,LONK
         IF(UEXT(I).GE.0.) THEN
           FM(I)=UEXT(I)*DTX
           ALF(I)=FM(I)/TM(I)
           TM(I)=TM(I)-FM(I)
         ENDIF
 141  CONTINUE
C
      DO 142 I=1,LONK
         ALFQ(I)=ALF(I)*ALF(I)
         ALF1(I)=1.-ALF(I)
         ALF1Q(I)=ALF1(I)*ALF1(I)
         ALF2(I)=ALF1(I)-ALF(I)
         ALF3(I)=ALF(I)*ALFQ(I)
         ALF4(I)=ALF1(I)*ALF1Q(I)
 142  CONTINUE
C
      DO 150 JV=1,NTRA
      DO 1500 I=1,LONK-1
C
         IF(UEXT(I).LT.0.) THEN
C
           F0 (I,JV)=ALF (I)* ( T0(I+1,JV)-ALF1(I)*
     +             ( TX(I+1,JV)-ALF2(I)*TXX(I+1,JV) ) )
           FX (I,JV)=ALFQ(I)*(TX(I+1,JV)-3.*ALF1(I)*TXX(I+1,JV))
           FXX(I,JV)=ALF3(I)*TXX(I+1,JV)
           FY (I,JV)=ALF (I)*(TY(I+1,JV)-ALF1(I)*TXY(I+1,JV))
           FZ (I,JV)=ALF (I)*(TZ(I+1,JV)-ALF1(I)*TXZ(I+1,JV))
           FXY(I,JV)=ALFQ(I)*TXY(I+1,JV)
           FXZ(I,JV)=ALFQ(I)*TXZ(I+1,JV)
           FYY(I,JV)=ALF (I)*TYY(I+1,JV)
           FYZ(I,JV)=ALF (I)*TYZ(I+1,JV)
           FZZ(I,JV)=ALF (I)*TZZ(I+1,JV)
C
           T0 (I+1,JV)=T0(I+1,JV)-F0(I,JV)
           TX (I+1,JV)=ALF1Q(I)*(TX(I+1,JV)+3.*ALF(I)*TXX(I+1,JV))
           TXX(I+1,JV)=ALF4(I)*TXX(I+1,JV)
           TY (I+1,JV)=TY (I+1,JV)-FY (I,JV)
           TZ (I+1,JV)=TZ (I+1,JV)-FZ (I,JV)
           TYY(I+1,JV)=TYY(I+1,JV)-FYY(I,JV)
           TYZ(I+1,JV)=TYZ(I+1,JV)-FYZ(I,JV)
           TZZ(I+1,JV)=TZZ(I+1,JV)-FZZ(I,JV)
           TXY(I+1,JV)=ALF1Q(I)*TXY(I+1,JV)
           TXZ(I+1,JV)=ALF1Q(I)*TXZ(I+1,JV)
C
         ENDIF
C
 1500 CONTINUE
 150  CONTINUE
C
      I=LONK
      IF(UEXT(I).LT.0.) THEN
C
        DO 151 JV=1,NTRA
C
           F0 (I,JV)=ALF (I)* ( T0(1,JV)-ALF1(I)*
     +             ( TX(1,JV)-ALF2(I)*TXX(1,JV) ) )
           FX (I,JV)=ALFQ(I)*(TX(1,JV)-3.*ALF1(I)*TXX(1,JV))
           FXX(I,JV)=ALF3(I)*TXX(1,JV)
           FY (I,JV)=ALF (I)*(TY(1,JV)-ALF1(I)*TXY(1,JV))
           FZ (I,JV)=ALF (I)*(TZ(1,JV)-ALF1(I)*TXZ(1,JV))
           FXY(I,JV)=ALFQ(I)*TXY(1,JV)
           FXZ(I,JV)=ALFQ(I)*TXZ(1,JV)
           FYY(I,JV)=ALF (I)*TYY(1,JV)
           FYZ(I,JV)=ALF (I)*TYZ(1,JV)
           FZZ(I,JV)=ALF (I)*TZZ(1,JV)
C
           T0 (1,JV)=T0(1,JV)-F0(I,JV)
           TX (1,JV)=ALF1Q(I)*(TX(1,JV)+3.*ALF(I)*TXX(1,JV))
           TXX(1,JV)=ALF4(I)*TXX(1,JV)
           TY (1,JV)=TY (1,JV)-FY (I,JV)
           TZ (1,JV)=TZ (1,JV)-FZ (I,JV)
           TYY(1,JV)=TYY(1,JV)-FYY(I,JV)
           TYZ(1,JV)=TYZ(1,JV)-FYZ(I,JV)
           TZZ(1,JV)=TZZ(1,JV)-FZZ(I,JV)
           TXY(1,JV)=ALF1Q(I)*TXY(1,JV)
           TXZ(1,JV)=ALF1Q(I)*TXZ(1,JV)
C
 151    CONTINUE
C
      ENDIF
C
      DO 152 JV=1,NTRA
      DO 1520 I=1,LONK
C
         IF(UEXT(I).GE.0.) THEN
C
           F0 (I,JV)=ALF (I)* ( T0(I,JV)+ALF1(I)*
     +             ( TX(I,JV)+ALF2(I)*TXX(I,JV) ) )
           FX (I,JV)=ALFQ(I)*(TX(I,JV)+3.*ALF1(I)*TXX(I,JV))
           FXX(I,JV)=ALF3(I)*TXX(I,JV)
           FY (I,JV)=ALF (I)*(TY(I,JV)+ALF1(I)*TXY(I,JV))
           FZ (I,JV)=ALF (I)*(TZ(I,JV)+ALF1(I)*TXZ(I,JV))
           FXY(I,JV)=ALFQ(I)*TXY(I,JV)
           FXZ(I,JV)=ALFQ(I)*TXZ(I,JV)
           FYY(I,JV)=ALF (I)*TYY(I,JV)
           FYZ(I,JV)=ALF (I)*TYZ(I,JV)
           FZZ(I,JV)=ALF (I)*TZZ(I,JV)
C
           T0 (I,JV)=T0(I,JV)-F0(I,JV)
           TX (I,JV)=ALF1Q(I)*(TX(I,JV)-3.*ALF(I)*TXX(I,JV))
           TXX(I,JV)=ALF4(I)*TXX(I,JV)
           TY (I,JV)=TY (I,JV)-FY (I,JV)
           TZ (I,JV)=TZ (I,JV)-FZ (I,JV)
           TYY(I,JV)=TYY(I,JV)-FYY(I,JV)
           TYZ(I,JV)=TYZ(I,JV)-FYZ(I,JV)
           TZZ(I,JV)=TZZ(I,JV)-FZZ(I,JV)
           TXY(I,JV)=ALF1Q(I)*TXY(I,JV)
           TXZ(I,JV)=ALF1Q(I)*TXZ(I,JV)
C
         ENDIF
C
 1520 CONTINUE
 152  CONTINUE
C
C  puts the temporary moments Fi into appropriate neighboring boxes
C
      DO 160 I=1,LONK
         IF(UEXT(I).LT.0.) THEN
           TM(I)=TM(I)+FM(I)
           ALF(I)=FM(I)/TM(I)
         ENDIF
 160  CONTINUE
C
      DO 161 I=1,LONK-1
         IF(UEXT(I).GE.0.) THEN
           TM(I+1)=TM(I+1)+FM(I)
           ALF(I)=FM(I)/TM(I+1)
         ENDIF
 161  CONTINUE
C
      I=LONK
      IF(UEXT(I).GE.0.) THEN
        TM(1)=TM(1)+FM(I)
        ALF(I)=FM(I)/TM(1)
      ENDIF
C
      DO 162 I=1,LONK
         ALF1(I)=1.-ALF(I)
         ALFQ(I)=ALF(I)*ALF(I)
         ALF1Q(I)=ALF1(I)*ALF1(I)
         ALF2(I)=ALF1(I)-ALF(I)
         ALF3(I)=ALF(I)*ALF1(I)
 162  CONTINUE
C
      DO 170 JV=1,NTRA
      DO 1700 I=1,LONK
C
         IF(UEXT(I).LT.0.) THEN
C
           TEMPTM=-ALF(I)*T0(I,JV)+ALF1(I)*F0(I,JV)
           T0 (I,JV)=T0(I,JV)+F0(I,JV)
           TXX(I,JV)=ALFQ(I)*FXX(I,JV)+ALF1Q(I)*TXX(I,JV)
     +          +5.*( ALF3(I)*(FX(I,JV)-TX(I,JV))+ALF2(I)*TEMPTM )
           TX (I,JV)=ALF (I)*FX (I,JV)+ALF1(I)*TX (I,JV)+3.*TEMPTM
           TXY(I,JV)=ALF (I)*FXY(I,JV)+ALF1(I)*TXY(I,JV)
     +          +3.*(ALF1(I)*FY (I,JV)-ALF (I)*TY (I,JV))
           TXZ(I,JV)=ALF (I)*FXZ(I,JV)+ALF1(I)*TXZ(I,JV)
     +          +3.*(ALF1(I)*FZ (I,JV)-ALF (I)*TZ (I,JV))
           TY (I,JV)=TY (I,JV)+FY (I,JV)
           TZ (I,JV)=TZ (I,JV)+FZ (I,JV)
           TYY(I,JV)=TYY(I,JV)+FYY(I,JV)
           TYZ(I,JV)=TYZ(I,JV)+FYZ(I,JV)
           TZZ(I,JV)=TZZ(I,JV)+FZZ(I,JV)
C
         ENDIF
C
 1700 CONTINUE
 170  CONTINUE
C
      DO 171 JV=1,NTRA
      DO 1710 I=1,LONK-1
C
         IF(UEXT(I).GE.0.) THEN
C
           TEMPTM=ALF(I)*T0(I+1,JV)-ALF1(I)*F0(I,JV)
           T0 (I+1,JV)=T0(I+1,JV)+F0(I,JV)
           TXX(I+1,JV)=ALFQ(I)*FXX(I,JV)+ALF1Q(I)*TXX(I+1,JV)
     +           +5.*( ALF3(I)*(TX(I+1,JV)-FX(I,JV))-ALF2(I)*TEMPTM )
           TX (I+1,JV)=ALF(I)*FX (I  ,JV)+ALF1(I)*TX (I+1,JV)+3.*TEMPTM
           TXY(I+1,JV)=ALF(I)*FXY(I  ,JV)+ALF1(I)*TXY(I+1,JV)
     +            +3.*(ALF(I)*TY (I+1,JV)-ALF1(I)*FY (I  ,JV))
           TXZ(I+1,JV)=ALF(I)*FXZ(I  ,JV)+ALF1(I)*TXZ(I+1,JV)
     +            +3.*(ALF(I)*TZ (I+1,JV)-ALF1(I)*FZ (I  ,JV))
           TY (I+1,JV)=TY (I+1,JV)+FY (I,JV)
           TZ (I+1,JV)=TZ (I+1,JV)+FZ (I,JV)
           TYY(I+1,JV)=TYY(I+1,JV)+FYY(I,JV)
           TYZ(I+1,JV)=TYZ(I+1,JV)+FYZ(I,JV)
           TZZ(I+1,JV)=TZZ(I+1,JV)+FZZ(I,JV)
C
         ENDIF
C
 1710 CONTINUE
 171  CONTINUE
C
      I=LONK
      IF(UEXT(I).GE.0.) THEN
        DO 172 JV=1,NTRA
           TEMPTM=ALF(I)*T0(1,JV)-ALF1(I)*F0(I,JV)
           T0 (1,JV)=T0(1,JV)+F0(I,JV)
           TXX(1,JV)=ALFQ(I)*FXX(I,JV)+ALF1Q(I)*TXX(1,JV)
     +         +5.*( ALF3(I)*(TX(1,JV)-FX(I,JV))-ALF2(I)*TEMPTM )
           TX (1,JV)=ALF(I)*FX(I,JV)+ALF1(I)*TX(1,JV)+3.*TEMPTM
           TXY(1,JV)=ALF(I)*FXY(I,JV)+ALF1(I)*TXY(1,JV)
     +          +3.*(ALF(I)*TY (1,JV)-ALF1(I)*FY (I,JV))
           TXZ(1,JV)=ALF(I)*FXZ(I,JV)+ALF1(I)*TXZ(1,JV)
     +          +3.*(ALF(I)*TZ (1,JV)-ALF1(I)*FZ (I,JV))
           TY (1,JV)=TY (1,JV)+FY (I,JV)
           TZ (1,JV)=TZ (1,JV)+FZ (I,JV)
           TYY(1,JV)=TYY(1,JV)+FYY(I,JV)
           TYZ(1,JV)=TYZ(1,JV)+FYZ(I,JV)
           TZZ(1,JV)=TZZ(1,JV)+FZZ(I,JV)
 172    CONTINUE
      ENDIF
C
C  retour aux mailles d'origine (passage des Tij aux Sij)
C
      IF(NUMK.GT.1) THEN
C
      DO 18 I2=1,NUMK
C
         DO 180 I=1,LONK
C
            I3=I2+(I-1)*NUMK
            SM(I3,K,L)=SMNEW(I3)
            ALF(I)=SMNEW(I3)/TM(I)
            TM(I)=TM(I)-SMNEW(I3)
C
            ALFQ(I)=ALF(I)*ALF(I)
            ALF1(I)=1.-ALF(I)
            ALF1Q(I)=ALF1(I)*ALF1(I)
            ALF2(I)=ALF1(I)-ALF(I)
            ALF3(I)=ALF(I)*ALFQ(I)
            ALF4(I)=ALF1(I)*ALF1Q(I)
C
 180     CONTINUE
C
         DO 181 JV=1,NTRA
         DO 181 I=1,LONK
C
            I3=I2+(I-1)*NUMK
            S0 (I3,K,L,JV)=ALF (I)* ( T0(I,JV)-ALF1(I)*
     +              ( TX(I,JV)-ALF2(I)*TXX(I,JV) ) )
            SSX (I3,K,L,JV)=ALFQ(I)*(TX(I,JV)-3.*ALF1(I)*TXX(I,JV))
            SSXX(I3,K,L,JV)=ALF3(I)*TXX(I,JV)
            SY (I3,K,L,JV)=ALF (I)*(TY(I,JV)-ALF1(I)*TXY(I,JV))
            SZ (I3,K,L,JV)=ALF (I)*(TZ(I,JV)-ALF1(I)*TXZ(I,JV))
            SSXY(I3,K,L,JV)=ALFQ(I)*TXY(I,JV)
            SSXZ(I3,K,L,JV)=ALFQ(I)*TXZ(I,JV)
            SYY(I3,K,L,JV)=ALF (I)*TYY(I,JV)
            SYZ(I3,K,L,JV)=ALF (I)*TYZ(I,JV)
            SZZ(I3,K,L,JV)=ALF (I)*TZZ(I,JV)
C
C   reajusts moments remaining in the box
C
            T0 (I,JV)=T0(I,JV)-S0(I3,K,L,JV)
            TX (I,JV)=ALF1Q(I)*(TX(I,JV)+3.*ALF(I)*TXX(I,JV))
            TXX(I,JV)=ALF4 (I)*TXX(I,JV)
            TY (I,JV)=TY (I,JV)-SY (I3,K,L,JV)
            TZ (I,JV)=TZ (I,JV)-SZ (I3,K,L,JV)
            TYY(I,JV)=TYY(I,JV)-SYY(I3,K,L,JV)
            TYZ(I,JV)=TYZ(I,JV)-SYZ(I3,K,L,JV)
            TZZ(I,JV)=TZZ(I,JV)-SZZ(I3,K,L,JV)
            TXY(I,JV)=ALF1Q(I)*TXY(I,JV)
            TXZ(I,JV)=ALF1Q(I)*TXZ(I,JV)
C
 181     CONTINUE
C
 18   CONTINUE
C
      ELSE
C
      DO 190 I=1,LON
         SM(I,K,L)=TM(I)
 190  CONTINUE
      DO 191 JV=1,NTRA
      DO 1910 I=1,LON
         S0 (I,K,L,JV)=T0 (I,JV)
         SSX (I,K,L,JV)=TX (I,JV)
         SY (I,K,L,JV)=TY (I,JV)
         SZ (I,K,L,JV)=TZ (I,JV)
         SSXX(I,K,L,JV)=TXX(I,JV)
         SSXY(I,K,L,JV)=TXY(I,JV)
         SSXZ(I,K,L,JV)=TXZ(I,JV)
         SYY(I,K,L,JV)=TYY(I,JV)
         SYZ(I,K,L,JV)=TYZ(I,JV)
         SZZ(I,K,L,JV)=TZZ(I,JV)
 1910 CONTINUE
 191  CONTINUE
C
      ENDIF
C
 1    CONTINUE
C
C ----------- AA Test en fin de ADVX ------ Controle des S*

c      DO 9999 l = 1, llm
c      DO 9999 j = 1, jjp1
c      DO 9999 i = 1, iip1
c	   IF (S0(i,j,l,ntra).lt.0..and.LIMIT) THEN
c           PRINT*, '-------------------'
c	        PRINT*, 'En fin de ADVXP'
c           PRINT*,'S0(',i,j,l,')=',S0(i,j,l,ntra)
c	        print*, 'SSX(',i,j,l,')=',SSX(i,j,l,ntra)
c           print*, 'SY(',i,j,l,')=',SY(i,j,l,ntra)
c       	print*, 'SZ(',i,j,l,')=',SZ(i,j,l,ntra)
c            WRITE (*,*) 'On arrete !! - pbl en fin de ADVXP'
c            STOP
c           ENDIF
c 9999 CONTINUE
c ---------- bouclage cyclique

      DO l = 1,llm
      DO j = 1,jjp1
         SM(iip1,j,l) = SM(1,j,l)
         S0(iip1,j,l,ntra) = S0(1,j,l,ntra)
     	 SSX(iip1,j,l,ntra) = SSX(1,j,l,ntra)
    	 SY(iip1,j,l,ntra) = SY(1,j,l,ntra)
    	 SZ(iip1,j,l,ntra) = SZ(1,j,l,ntra)
      END DO
      END DO

C ----------- qqtite totale de traceur dans tte l'atmosphere
      DO l = 1, llm
      DO j = 1, jjp1
      DO i = 1, iim
        sqf = sqf + S0(i,j,l,ntra)
      END DO
      END DO
      END DO

      PRINT*,'------ DIAG DANS ADVX2 - SORTIE -----'
      PRINT*,'sqf=',sqf
c-------------------------------------------------------------
      RETURN
      END
