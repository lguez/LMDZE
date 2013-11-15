!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/advx.F,v 1.2 2005/05/25 13:10:09 fairhead Exp $
!
      SUBROUTINE  advx(limit,dtx,pbaru,sm,s0,
     $     sx,sy,sz,lati,latf)
      use dimens_m
      use paramet_m
      use comconst
      use disvert_m
      IMPLICIT NONE

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                C
C  first-order moments (FOM) advection of tracer in X direction  C
C                                                                C
C  Source : Pascal Simon (Meteo,CNRM)                            C
C  Adaptation : A.Armengaud (LGGE) juin 94                       C
C                                                                C
C  limit,dtx,pbaru,pbarv,sm,s0,sx,sy,sz                       C
C  sont des arguments d'entree pour le s-pg...                   C
C                                                                C
C  sm,s0,sx,sy,sz                                                C
C  sont les arguments de sortie pour le s-pg                     C
C								 C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  parametres principaux du modele
C

C  Arguments :
C  -----------
C  dtx : frequence fictive d'appel du transport 
C  pbaru, pbarv : flux de masse en x et y en Pa.m2.s-1

       INTEGER ntra
       PARAMETER (ntra = 1)

C ATTENTION partout ou on trouve ntra, insertion de boucle
C           possible dans l'avenir.

      REAL dtx
      REAL, intent(in):: pbaru ( iip1,jjp1,llm )

C  moments: SM  total mass in each grid box
C           S0  mass of tracer in each grid box
C           Si  1rst order moment in i direction
C
      REAL SM(iip1,jjp1,llm),S0(iip1,jjp1,llm,ntra)
      REAL sx(iip1,jjp1,llm,ntra)
     $    ,sy(iip1,jjp1,llm,ntra)
      REAL sz(iip1,jjp1,llm,ntra)

C  Local :
C  ------- 

C  mass fluxes across the boundaries (UGRI,VGRI,WGRI)
C  mass fluxes in kg
C  declaration :

      REAL UGRI(iip1,jjp1,llm)

C  Rem : VGRI et WGRI ne sont pas utilises dans 
C  cette subroutine ( advection en x uniquement )
C
C  Ti are the moments for the current latitude and level
C
      REAL TM(iim)
      REAL T0(iim,ntra),TX(iim,ntra)
      REAL TY(iim,ntra),TZ(iim,ntra)
      REAL TEMPTM                ! just a temporary variable
C
C  the moments F are similarly defined and used as temporary
C  storage for portions of the grid boxes in transit
C
      REAL FM(iim)
      REAL F0(iim,ntra),FX(iim,ntra)
      REAL FY(iim,ntra),FZ(iim,ntra)
C
C  work arrays
C
      REAL ALF(iim),ALF1(iim),ALFQ(iim),ALF1Q(iim)
C
      REAL SMNEW(iim),UEXT(iim)
C
      REAL sqi,sqf

      LOGICAL LIMIT
      INTEGER NUM(jjp1),LONK,NUMK
      INTEGER lon,lati,latf,niv
      INTEGER i,i2,i3,j,jv,l,k,itrac 

      lon = iim 
      niv = llm 

C *** Test de passage d'arguments ******


C  -------------------------------------
      DO 300 j = 1,jjp1 
         NUM(j) = 1
  300 CONTINUE
      sqi = 0.
      sqf = 0.

      DO l = 1,llm
         DO j = 1,jjp1
            DO i = 1,iim
cIM 240305            sqi = sqi + S0(i,j,l,9)
               sqi = sqi + S0(i,j,l,ntra)
            ENDDO
         ENDDO
      ENDDO
      PRINT*,'-------- DIAG DANS ADVX - ENTREE ---------'
      PRINT*,'sqi=',sqi


C  Interface : adaptation nouveau modele
C  -------------------------------------
C
C  ---------------------------------------------------------
C  Conversion des flux de masses en kg/s
C  pbaru est en N/s d'ou :
C  ugri est en kg/s

      DO 500 l = 1,llm
         DO 500 j = 1,jjm+1
            DO 500 i = 1,iip1  
C            ugri (i,j,llm+1-l) = pbaru (i,j,l) * ( dsig(l) / g )
             ugri (i,j,llm+1-l) = pbaru (i,j,l)
  500 CONTINUE


C  ---------------------------------------------------------
C  ---------------------------------------------------------
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
         T0(I,JV)=0.
         TX(I,JV)=0.
         TY(I,JV)=0.
         TZ(I,JV)=0.
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
 113     CONTINUE
C
         DO  JV=1,NTRA
         DO  I=1,LONK
            I3=(I-1)*NUMK+I2
            TEMPTM=-ALF(I)*T0(I,JV)+ALF1(I)
     $          *S0(I3,K,L,JV)
            T0(I,JV)=T0(I,JV)+S0(I3,K,L,JV)
            TX(I,JV)=ALF(I)  *sx(I3,K,L,JV)+
     $       ALF1(I)*TX(I,JV) +3.*TEMPTM
            TY(I,JV)=TY(I,JV)+sy(I3,K,L,JV)
            TZ(I,JV)=TZ(I,JV)+sz(I3,K,L,JV)
         ENDDO 
         ENDDO
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
         T0(I,JV)=S0(I,K,L,JV)
         TX(I,JV)=sx(I,K,L,JV)
         TY(I,JV)=sy(I,K,L,JV)
         TZ(I,JV)=sz(I,K,L,JV)
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
        TX(I,JV)=SIGN(AMIN1(AMAX1(T0(I,JV),0.),ABS(TX(I,JV))),TX(I,JV))
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
 142  CONTINUE
C
      DO 150 JV=1,NTRA
      DO 1500 I=1,LONK-1
C
         IF(UEXT(I).LT.0.) THEN
C
           F0(I,JV)=ALF (I)* ( T0(I+1,JV)-ALF1(I)*TX(I+1,JV) )
           FX(I,JV)=ALFQ(I)*TX(I+1,JV)
           FY(I,JV)=ALF (I)*TY(I+1,JV)
           FZ(I,JV)=ALF (I)*TZ(I+1,JV)
C
           T0(I+1,JV)=T0(I+1,JV)-F0(I,JV)
           TX(I+1,JV)=ALF1Q(I)*TX(I+1,JV)
           TY(I+1,JV)=TY(I+1,JV)-FY(I,JV)
           TZ(I+1,JV)=TZ(I+1,JV)-FZ(I,JV)
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
           F0 (I,JV)=ALF (I)* ( T0(1,JV)-ALF1(I)*TX(1,JV) )
           FX (I,JV)=ALFQ(I)*TX(1,JV)
           FY (I,JV)=ALF (I)*TY(1,JV)
           FZ (I,JV)=ALF (I)*TZ(1,JV)
C
           T0(1,JV)=T0(1,JV)-F0(I,JV)
           TX(1,JV)=ALF1Q(I)*TX(1,JV)
           TY(1,JV)=TY(1,JV)-FY(I,JV)
           TZ(1,JV)=TZ(1,JV)-FZ(I,JV)
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
           F0(I,JV)=ALF (I)* ( T0(I,JV)+ALF1(I)*TX(I,JV) )
           FX(I,JV)=ALFQ(I)*TX(I,JV)
           FY(I,JV)=ALF (I)*TY(I,JV)
           FZ(I,JV)=ALF (I)*TZ(I,JV)
C
           T0(I,JV)=T0(I,JV)-F0(I,JV)
           TX(I,JV)=ALF1Q(I)*TX(I,JV)
           TY(I,JV)=TY(I,JV)-FY(I,JV)
           TZ(I,JV)=TZ(I,JV)-FZ(I,JV)
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
 162  CONTINUE
C
      DO 170 JV=1,NTRA
      DO 1700 I=1,LONK
C
         IF(UEXT(I).LT.0.) THEN
C
           TEMPTM=-ALF(I)*T0(I,JV)+ALF1(I)*F0(I,JV)
           T0(I,JV)=T0(I,JV)+F0(I,JV)
           TX(I,JV)=ALF(I)*FX(I,JV)+ALF1(I)*TX(I,JV)+3.*TEMPTM
           TY(I,JV)=TY(I,JV)+FY(I,JV)
           TZ(I,JV)=TZ(I,JV)+FZ(I,JV)
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
           T0(I+1,JV)=T0(I+1,JV)+F0(I,JV)
           TX(I+1,JV)=ALF(I)*FX(I,JV)+ALF1(I)*TX(I+1,JV)+3.*TEMPTM
           TY(I+1,JV)=TY(I+1,JV)+FY(I,JV)
           TZ(I+1,JV)=TZ(I+1,JV)+FZ(I,JV)
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
           T0(1,JV)=T0(1,JV)+F0(I,JV)
           TX(1,JV)=ALF(I)*FX(I,JV)+ALF1(I)*TX(1,JV)+3.*TEMPTM
           TY(1,JV)=TY(1,JV)+FY(I,JV)
           TZ(1,JV)=TZ(1,JV)+FZ(I,JV)
 172    CONTINUE
      ENDIF
C
C  retour aux mailles d'origine (passage des Tij aux Sij)
C
      IF(NUMK.GT.1) THEN
C
      DO 180 I2=1,NUMK
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
C
 180     CONTINUE
C
         DO  JV=1,NTRA
         DO  I=1,LONK
C
            I3=I2+(I-1)*NUMK
            S0(I3,K,L,JV)=ALF (I)
     $       * (T0(I,JV)-ALF1(I)*TX(I,JV))
            sx(I3,K,L,JV)=ALFQ(I)*TX(I,JV)
            sy(I3,K,L,JV)=ALF (I)*TY(I,JV)
            sz(I3,K,L,JV)=ALF (I)*TZ(I,JV)
C
C   reajusts moments remaining in the box
C
            T0(I,JV)=T0(I,JV)-S0(I3,K,L,JV)
            TX(I,JV)=ALF1Q(I)*TX(I,JV)
            TY(I,JV)=TY(I,JV)-sy(I3,K,L,JV)
            TZ(I,JV)=TZ(I,JV)-sz(I3,K,L,JV)
          ENDDO
          ENDDO
C
C
      ELSE
C
      DO 190 I=1,LON
         SM(I,K,L)=TM(I)
 190  CONTINUE
      DO 191 JV=1,NTRA
      DO 1910 I=1,LON
         S0(I,K,L,JV)=T0(I,JV)
         sx(I,K,L,JV)=TX(I,JV)
         sy(I,K,L,JV)=TY(I,JV)
         sz(I,K,L,JV)=TZ(I,JV)
 1910 CONTINUE
 191  CONTINUE
C
      ENDIF
C
 1    CONTINUE
C
C ---------- bouclage cyclique 
      DO itrac=1,ntra
      DO l = 1,llm
        DO j = lati,latf
           SM(iip1,j,l) = SM(1,j,l)
           S0(iip1,j,l,itrac) = S0(1,j,l,itrac)
           sx(iip1,j,l,itrac) = sx(1,j,l,itrac)
           sy(iip1,j,l,itrac) = sy(1,j,l,itrac)
           sz(iip1,j,l,itrac) = sz(1,j,l,itrac)
        END DO
      END DO
      ENDDO 

c ----------- qqtite totale de traceur dans tte l'atmosphere
      DO l = 1, llm
        DO j = 1, jjp1
          DO i = 1, iim
cIM 240405          sqf = sqf + S0(i,j,l,9)
             sqf = sqf + S0(i,j,l,ntra)
          END DO  
        END DO
      END DO
c
      PRINT*,'------ DIAG DANS ADVX - SORTIE -----'
      PRINT*,'sqf=',sqf
c-------------

      RETURN
      END
C_________________________________________________________________
C_________________________________________________________________
