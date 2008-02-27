!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/advyp.F,v 1.1.1.1 2004/05/19 12:53:06 lmdzadmin Exp $
!
      SUBROUTINE ADVYP(LIMIT,DTY,PBARV,SM,S0,SSX,SY,SZ
     .                 ,SSXX,SSXY,SSXZ,SYY,SYZ,SZZ,ntra )
      use dimens_m
      use comconst
      use paramet_m
      use comvert
      use comgeom
      IMPLICIT NONE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                 C
C  second-order moments (SOM) advection of tracer in Y direction  C
C                                                                 C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                C
C  Source : Pascal Simon ( Meteo, CNRM )			 C
C  Adaptation : A.A. (LGGE) 					 C
C  Derniere Modif : 19/10/95 LAST
C								 C
C  sont les arguments d'entree pour le s-pg			 C
C								 C
C  argument de sortie du s-pg					 C
C								 C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Rem : Probleme aux poles il faut reecrire ce cas specifique
C        Attention au sens de l'indexation 
C
C  parametres principaux du modele
C
C
 
C  Arguments :
C  ----------
C  dty : frequence fictive d'appel du transport
C  parbu,pbarv : flux de masse en x et y en Pa.m2.s-1

      INTEGER lon,lat,niv
      INTEGER i,j,jv,k,kp,l
      INTEGER ntra
C      PARAMETER (ntra = 1)

      REAL dty
      REAL pbarv ( iip1,jjm, llm )

C  moments: SM  total mass in each grid box
C           S0  mass of tracer in each grid box
C           Si  1rst order moment in i direction
C
      REAL SM(iip1,jjp1,llm)
     +    ,S0(iip1,jjp1,llm,ntra)
      REAL SSX(iip1,jjp1,llm,ntra)
     +    ,SY(iip1,jjp1,llm,ntra)
     +    ,SZ(iip1,jjp1,llm,ntra)
     +    ,SSXX(iip1,jjp1,llm,ntra)
     +    ,SSXY(iip1,jjp1,llm,ntra)
     +    ,SSXZ(iip1,jjp1,llm,ntra)
     +    ,SYY(iip1,jjp1,llm,ntra)
     +    ,SYZ(iip1,jjp1,llm,ntra)
     +    ,SZZ(iip1,jjp1,llm,ntra)
C
C  Local :
C  ------- 

C  mass fluxes across the boundaries (UGRI,VGRI,WGRI)
C  mass fluxes in kg
C  declaration :

      REAL VGRI(iip1,0:jjp1,llm)

C  Rem : UGRI et WGRI ne sont pas utilises dans 
C  cette subroutine ( advection en y uniquement )
C  Rem 2 :le dimensionnement de VGRI depend de celui de pbarv
C
C  the moments F are similarly defined and used as temporary
C  storage for portions of the grid boxes in transit
C
C  the moments Fij are used as temporary storage for
C  portions of the grid boxes in transit at the current level
C
C  work arrays
C
C
      REAL F0(iim,0:jjp1,ntra),FM(iim,0:jjp1)
      REAL FX(iim,jjm,ntra),FY(iim,jjm,ntra)
      REAL FZ(iim,jjm,ntra)
      REAL FXX(iim,jjm,ntra),FXY(iim,jjm,ntra)
      REAL FXZ(iim,jjm,ntra),FYY(iim,jjm,ntra)
      REAL FYZ(iim,jjm,ntra),FZZ(iim,jjm,ntra)
      REAL S00(ntra)
      REAL SM0             ! Just temporal variable
C
C  work arrays
C
      REAL ALF(iim,0:jjp1),ALF1(iim,0:jjp1)
      REAL ALFQ(iim,0:jjp1),ALF1Q(iim,0:jjp1)
      REAL ALF2(iim,0:jjp1),ALF3(iim,0:jjp1)
      REAL ALF4(iim,0:jjp1)
      REAL TEMPTM          ! Just temporal variable
      REAL SLPMAX,S1MAX,S1NEW,S2NEW
c
C  Special pour poles 
c
      REAL sbms,sfms,sfzs,sbmn,sfmn,sfzn
      REAL sns0(ntra),snsz(ntra),snsm
      REAL qy1(iim,llm,ntra),qylat(iim,llm,ntra)
      REAL cx1(llm,ntra), cxLAT(llm,ntra)
      REAL cy1(llm,ntra), cyLAT(llm,ntra)
      REAL z1(iim), zcos(iim), zsin(iim)
      REAL SSUM
      EXTERNAL SSUM
C
      REAL sqi,sqf
      LOGICAL LIMIT

      lon = iim         ! rem : Il est possible qu'un pbl. arrive ici
      lat = jjp1        ! a cause des dim. differentes entre les
      niv = llm         !       tab. S et VGRI 
                    
c-----------------------------------------------------------------
C initialisations

      sbms = 0.
      sfms = 0.
      sfzs = 0.
      sbmn = 0.
      sfmn = 0.
      sfzn = 0.

c-----------------------------------------------------------------
C *** Test : diag de la qtite totale de traceur dans
C            l'atmosphere avant l'advection en Y
c 
      sqi = 0.
      sqf = 0.

      DO l = 1,llm
         DO j = 1,jjp1
           DO i = 1,iim
              sqi = sqi + S0(i,j,l,ntra)
           END DO
         END DO
      END DO
      PRINT*,'---------- DIAG DANS ADVY - ENTREE --------'
      PRINT*,'sqi=',sqi

c-----------------------------------------------------------------
C  Interface : adaptation nouveau modele
C  -------------------------------------
C
C  Conversion des flux de masses en kg
C-AA 20/10/94  le signe -1 est necessaire car indexation opposee

      DO 500 l = 1,llm
         DO 500 j = 1,jjm
            DO 500 i = 1,iip1  
            vgri (i,j,llm+1-l)=-1.*pbarv (i,j,l)
  500 CONTINUE

CAA Initialisation de flux fictifs aux bords sup. des boites pol.

      DO l = 1,llm
         DO i = 1,iip1  
             vgri(i,0,l) = 0.
             vgri(i,jjp1,l) = 0.
         ENDDO
      ENDDO
c
c----------------- START HERE -----------------------
C  boucle sur les niveaux
C
      DO 1 L=1,NIV
C
C  place limits on appropriate moments before transport
C      (if flux-limiting is to be applied)
C
      IF(.NOT.LIMIT) GO TO 11
C
      DO 10 JV=1,NTRA
      DO 10 K=1,LAT
      DO 100 I=1,LON
         IF(S0(I,K,L,JV).GT.0.) THEN
           SLPMAX=AMAX1(S0(I,K,L,JV),0.)
           S1MAX=1.5*SLPMAX
           S1NEW=AMIN1(S1MAX,AMAX1(-S1MAX,SY(I,K,L,JV)))
           S2NEW=AMIN1( 2.*SLPMAX-ABS(S1NEW)/3. ,
     +                  AMAX1(ABS(S1NEW)-SLPMAX,SYY(I,K,L,JV)) )
           SY (I,K,L,JV)=S1NEW
           SYY(I,K,L,JV)=S2NEW
       SSXY(I,K,L,JV)=AMIN1(SLPMAX,AMAX1(-SLPMAX,SSXY(I,K,L,JV)))
       SYZ(I,K,L,JV)=AMIN1(SLPMAX,AMAX1(-SLPMAX,SYZ(I,K,L,JV)))
         ELSE
           SY (I,K,L,JV)=0.
           SYY(I,K,L,JV)=0.
           SSXY(I,K,L,JV)=0.
           SYZ(I,K,L,JV)=0.
         ENDIF
 100  CONTINUE
 10   CONTINUE
C
 11   CONTINUE
C
C  le flux a travers le pole Nord est traite separement
C
      SM0=0.
      DO 20 JV=1,NTRA
         S00(JV)=0.
 20   CONTINUE
C
      DO 21 I=1,LON
C
         IF(VGRI(I,0,L).LE.0.) THEN
           FM(I,0)=-VGRI(I,0,L)*DTY
           ALF(I,0)=FM(I,0)/SM(I,1,L)
           SM(I,1,L)=SM(I,1,L)-FM(I,0)
           SM0=SM0+FM(I,0)
         ENDIF
C
         ALFQ(I,0)=ALF(I,0)*ALF(I,0)
         ALF1(I,0)=1.-ALF(I,0)
         ALF1Q(I,0)=ALF1(I,0)*ALF1(I,0)
         ALF2(I,0)=ALF1(I,0)-ALF(I,0)
         ALF3(I,0)=ALF(I,0)*ALFQ(I,0)
         ALF4(I,0)=ALF1(I,0)*ALF1Q(I,0)
C
 21   CONTINUE
c     print*,'ADVYP 21'
C
      DO 22 JV=1,NTRA
      DO 220 I=1,LON
C
         IF(VGRI(I,0,L).LE.0.) THEN
C
           F0(I,0,JV)=ALF(I,0)* ( S0(I,1,L,JV)-ALF1(I,0)*
     +        ( SY(I,1,L,JV)-ALF2(I,0)*SYY(I,1,L,JV) ) )
C
           S00(JV)=S00(JV)+F0(I,0,JV)
           S0 (I,1,L,JV)=S0(I,1,L,JV)-F0(I,0,JV)
           SY (I,1,L,JV)=ALF1Q(I,0)*
     +            (SY(I,1,L,JV)+3.*ALF(I,0)*SYY(I,1,L,JV))
           SYY(I,1,L,JV)=ALF4 (I,0)*SYY(I,1,L,JV)
           SSX (I,1,L,JV)=ALF1 (I,0)*
     +            (SSX(I,1,L,JV)+ALF(I,0)*SSXY(I,1,L,JV) )
           SZ (I,1,L,JV)=ALF1 (I,0)*
     +            (SZ(I,1,L,JV)+ALF(I,0)*SSXZ(I,1,L,JV) )
           SSXX(I,1,L,JV)=ALF1 (I,0)*SSXX(I,1,L,JV)
           SSXZ(I,1,L,JV)=ALF1 (I,0)*SSXZ(I,1,L,JV)
           SZZ(I,1,L,JV)=ALF1 (I,0)*SZZ(I,1,L,JV)
           SSXY(I,1,L,JV)=ALF1Q(I,0)*SSXY(I,1,L,JV)
           SYZ(I,1,L,JV)=ALF1Q(I,0)*SYZ(I,1,L,JV)
C
         ENDIF
C
 220  CONTINUE
 22   CONTINUE
C
      DO 23 I=1,LON
         IF(VGRI(I,0,L).GT.0.) THEN
           FM(I,0)=VGRI(I,0,L)*DTY
           ALF(I,0)=FM(I,0)/SM0
         ENDIF
 23   CONTINUE
C
      DO 24 JV=1,NTRA
      DO 240 I=1,LON
         IF(VGRI(I,0,L).GT.0.) THEN
           F0(I,0,JV)=ALF(I,0)*S00(JV)
         ENDIF
 240  CONTINUE
 24   CONTINUE
C
C  puts the temporary moments Fi into appropriate neighboring boxes
C
c     print*,'av ADVYP 25'
      DO 25 I=1,LON
C
         IF(VGRI(I,0,L).GT.0.) THEN
           SM(I,1,L)=SM(I,1,L)+FM(I,0)
           ALF(I,0)=FM(I,0)/SM(I,1,L)
         ENDIF
C
         ALFQ(I,0)=ALF(I,0)*ALF(I,0)
         ALF1(I,0)=1.-ALF(I,0)
         ALF1Q(I,0)=ALF1(I,0)*ALF1(I,0)
         ALF2(I,0)=ALF1(I,0)-ALF(I,0)
         ALF3(I,0)=ALF1(I,0)*ALF(I,0)
C
 25   CONTINUE
c     print*,'av ADVYP 25'
C
      DO 26 JV=1,NTRA
      DO 260 I=1,LON
C
         IF(VGRI(I,0,L).GT.0.) THEN
C
         TEMPTM=ALF(I,0)*S0(I,1,L,JV)-ALF1(I,0)*F0(I,0,JV)
         S0 (I,1,L,JV)=S0(I,1,L,JV)+F0(I,0,JV)
         SYY(I,1,L,JV)=ALF1Q(I,0)*SYY(I,1,L,JV)
     +        +5.*( ALF3 (I,0)*SY (I,1,L,JV)-ALF2(I,0)*TEMPTM )
         SY (I,1,L,JV)=ALF1 (I,0)*SY (I,1,L,JV)+3.*TEMPTM
      SSXY(I,1,L,JV)=ALF1 (I,0)*SSXY(I,1,L,JV)+3.*ALF(I,0)*SSX(I,1,L,JV)
      SYZ(I,1,L,JV)=ALF1 (I,0)*SYZ(I,1,L,JV)+3.*ALF(I,0)*SZ(I,1,L,JV)
C
         ENDIF
C
 260  CONTINUE
 26   CONTINUE
C
C  calculate flux and moments between adjacent boxes
C  1- create temporary moments/masses for partial boxes in transit
C  2- reajusts moments remaining in the box
C
C  flux from KP to K if V(K).lt.0 and from K to KP if V(K).gt.0
C
c     print*,'av ADVYP 30'
      DO 30 K=1,LAT-1
      KP=K+1
      DO 300 I=1,LON
C
         IF(VGRI(I,K,L).LT.0.) THEN
           FM(I,K)=-VGRI(I,K,L)*DTY
           ALF(I,K)=FM(I,K)/SM(I,KP,L)
           SM(I,KP,L)=SM(I,KP,L)-FM(I,K)
         ELSE
           FM(I,K)=VGRI(I,K,L)*DTY
           ALF(I,K)=FM(I,K)/SM(I,K,L)
           SM(I,K,L)=SM(I,K,L)-FM(I,K)
         ENDIF
C
         ALFQ(I,K)=ALF(I,K)*ALF(I,K)
         ALF1(I,K)=1.-ALF(I,K)
         ALF1Q(I,K)=ALF1(I,K)*ALF1(I,K)
         ALF2(I,K)=ALF1(I,K)-ALF(I,K)
         ALF3(I,K)=ALF(I,K)*ALFQ(I,K)
         ALF4(I,K)=ALF1(I,K)*ALF1Q(I,K)
C
 300  CONTINUE
 30   CONTINUE
c     print*,'ap ADVYP 30'
C
      DO 31 JV=1,NTRA
      DO 31 K=1,LAT-1
      KP=K+1
      DO 310 I=1,LON
C
         IF(VGRI(I,K,L).LT.0.) THEN
C
           F0 (I,K,JV)=ALF (I,K)* ( S0(I,KP,L,JV)-ALF1(I,K)*
     +        ( SY(I,KP,L,JV)-ALF2(I,K)*SYY(I,KP,L,JV) ) )
           FY (I,K,JV)=ALFQ(I,K)*
     +                 (SY(I,KP,L,JV)-3.*ALF1(I,K)*SYY(I,KP,L,JV))
           FYY(I,K,JV)=ALF3(I,K)*SYY(I,KP,L,JV)
           FX (I,K,JV)=ALF (I,K)*
     +                 (SSX(I,KP,L,JV)-ALF1(I,K)*SSXY(I,KP,L,JV))
           FZ (I,K,JV)=ALF (I,K)*
     +                 (SZ(I,KP,L,JV)-ALF1(I,K)*SYZ(I,KP,L,JV))
           FXY(I,K,JV)=ALFQ(I,K)*SSXY(I,KP,L,JV)
           FYZ(I,K,JV)=ALFQ(I,K)*SYZ(I,KP,L,JV)
           FXX(I,K,JV)=ALF (I,K)*SSXX(I,KP,L,JV)
           FXZ(I,K,JV)=ALF (I,K)*SSXZ(I,KP,L,JV)
           FZZ(I,K,JV)=ALF (I,K)*SZZ(I,KP,L,JV)
C
           S0 (I,KP,L,JV)=S0(I,KP,L,JV)-F0(I,K,JV)
           SY (I,KP,L,JV)=ALF1Q(I,K)*
     +                 (SY(I,KP,L,JV)+3.*ALF(I,K)*SYY(I,KP,L,JV))
           SYY(I,KP,L,JV)=ALF4(I,K)*SYY(I,KP,L,JV)
           SSX (I,KP,L,JV)=SSX (I,KP,L,JV)-FX (I,K,JV)
           SZ (I,KP,L,JV)=SZ (I,KP,L,JV)-FZ (I,K,JV)
           SSXX(I,KP,L,JV)=SSXX(I,KP,L,JV)-FXX(I,K,JV)
           SSXZ(I,KP,L,JV)=SSXZ(I,KP,L,JV)-FXZ(I,K,JV)
           SZZ(I,KP,L,JV)=SZZ(I,KP,L,JV)-FZZ(I,K,JV)
           SSXY(I,KP,L,JV)=ALF1Q(I,K)*SSXY(I,KP,L,JV)
           SYZ(I,KP,L,JV)=ALF1Q(I,K)*SYZ(I,KP,L,JV)
C
         ELSE
C
           F0 (I,K,JV)=ALF (I,K)* ( S0(I,K,L,JV)+ALF1(I,K)*
     +             ( SY(I,K,L,JV)+ALF2(I,K)*SYY(I,K,L,JV) ) )
           FY (I,K,JV)=ALFQ(I,K)*
     +                 (SY(I,K,L,JV)+3.*ALF1(I,K)*SYY(I,K,L,JV))
           FYY(I,K,JV)=ALF3(I,K)*SYY(I,K,L,JV)
      FX (I,K,JV)=ALF (I,K)*(SSX(I,K,L,JV)+ALF1(I,K)*SSXY(I,K,L,JV))
      FZ (I,K,JV)=ALF (I,K)*(SZ(I,K,L,JV)+ALF1(I,K)*SYZ(I,K,L,JV))
           FXY(I,K,JV)=ALFQ(I,K)*SSXY(I,K,L,JV)
           FYZ(I,K,JV)=ALFQ(I,K)*SYZ(I,K,L,JV)
           FXX(I,K,JV)=ALF (I,K)*SSXX(I,K,L,JV)
           FXZ(I,K,JV)=ALF (I,K)*SSXZ(I,K,L,JV)
           FZZ(I,K,JV)=ALF (I,K)*SZZ(I,K,L,JV)
C
           S0 (I,K,L,JV)=S0 (I,K,L,JV)-F0 (I,K,JV)
           SY (I,K,L,JV)=ALF1Q(I,K)*
     +                  (SY(I,K,L,JV)-3.*ALF(I,K)*SYY(I,K,L,JV))
           SYY(I,K,L,JV)=ALF4(I,K)*SYY(I,K,L,JV)
           SSX (I,K,L,JV)=SSX (I,K,L,JV)-FX (I,K,JV)
           SZ (I,K,L,JV)=SZ (I,K,L,JV)-FZ (I,K,JV)
           SSXX(I,K,L,JV)=SSXX(I,K,L,JV)-FXX(I,K,JV)
           SSXZ(I,K,L,JV)=SSXZ(I,K,L,JV)-FXZ(I,K,JV)
           SZZ(I,K,L,JV)=SZZ(I,K,L,JV)-FZZ(I,K,JV)
           SSXY(I,K,L,JV)=ALF1Q(I,K)*SSXY(I,K,L,JV)
           SYZ(I,K,L,JV)=ALF1Q(I,K)*SYZ(I,K,L,JV)
C
         ENDIF
C
 310  CONTINUE
 31   CONTINUE
c     print*,'ap ADVYP 31'
C
C  puts the temporary moments Fi into appropriate neighboring boxes
C
      DO 32 K=1,LAT-1
      KP=K+1
      DO 320 I=1,LON
C
         IF(VGRI(I,K,L).LT.0.) THEN
           SM(I,K,L)=SM(I,K,L)+FM(I,K)
           ALF(I,K)=FM(I,K)/SM(I,K,L)
         ELSE
           SM(I,KP,L)=SM(I,KP,L)+FM(I,K)
           ALF(I,K)=FM(I,K)/SM(I,KP,L)
         ENDIF
C
         ALFQ(I,K)=ALF(I,K)*ALF(I,K)
         ALF1(I,K)=1.-ALF(I,K)
         ALF1Q(I,K)=ALF1(I,K)*ALF1(I,K)
         ALF2(I,K)=ALF1(I,K)-ALF(I,K)
         ALF3(I,K)=ALF1(I,K)*ALF(I,K)
C
 320  CONTINUE
 32   CONTINUE
c     print*,'ap ADVYP 32'
C
      DO 33 JV=1,NTRA
      DO 33 K=1,LAT-1
      KP=K+1
      DO 330 I=1,LON
C
         IF(VGRI(I,K,L).LT.0.) THEN
C
         TEMPTM=-ALF(I,K)*S0(I,K,L,JV)+ALF1(I,K)*F0(I,K,JV)
         S0 (I,K,L,JV)=S0(I,K,L,JV)+F0(I,K,JV)
       SYY(I,K,L,JV)=ALFQ(I,K)*FYY(I,K,JV)+ALF1Q(I,K)*SYY(I,K,L,JV)
     +  +5.*( ALF3(I,K)*(FY(I,K,JV)-SY(I,K,L,JV))+ALF2(I,K)*TEMPTM )
         SY (I,K,L,JV)=ALF(I,K)*FY(I,K,JV)+ALF1(I,K)*SY(I,K,L,JV)
     +            +3.*TEMPTM
       SSXY(I,K,L,JV)=ALF (I,K)*FXY(I,K,JV)+ALF1(I,K)*SSXY(I,K,L,JV)
     +         +3.*(ALF1(I,K)*FX (I,K,JV)-ALF (I,K)*SSX (I,K,L,JV))
       SYZ(I,K,L,JV)=ALF (I,K)*FYZ(I,K,JV)+ALF1(I,K)*SYZ(I,K,L,JV)
     +         +3.*(ALF1(I,K)*FZ (I,K,JV)-ALF (I,K)*SZ (I,K,L,JV))
         SSX (I,K,L,JV)=SSX (I,K,L,JV)+FX (I,K,JV)
         SZ (I,K,L,JV)=SZ (I,K,L,JV)+FZ (I,K,JV)
         SSXX(I,K,L,JV)=SSXX(I,K,L,JV)+FXX(I,K,JV)
         SSXZ(I,K,L,JV)=SSXZ(I,K,L,JV)+FXZ(I,K,JV)
         SZZ(I,K,L,JV)=SZZ(I,K,L,JV)+FZZ(I,K,JV)
C
         ELSE
C
         TEMPTM=ALF(I,K)*S0(I,KP,L,JV)-ALF1(I,K)*F0(I,K,JV)
         S0 (I,KP,L,JV)=S0(I,KP,L,JV)+F0(I,K,JV)
       SYY(I,KP,L,JV)=ALFQ(I,K)*FYY(I,K,JV)+ALF1Q(I,K)*SYY(I,KP,L,JV)
     +  +5.*( ALF3(I,K)*(SY(I,KP,L,JV)-FY(I,K,JV))-ALF2(I,K)*TEMPTM )
         SY (I,KP,L,JV)=ALF(I,K)*FY(I,K,JV)+ALF1(I,K)*SY(I,KP,L,JV)
     +                 +3.*TEMPTM
       SSXY(I,KP,L,JV)=ALF(I,K)*FXY(I,K,JV)+ALF1(I,K)*SSXY(I,KP,L,JV)
     +             +3.*(ALF(I,K)*SSX(I,KP,L,JV)-ALF1(I,K)*FX(I,K,JV))
         SYZ(I,KP,L,JV)=ALF(I,K)*FYZ(I,K,JV)+ALF1(I,K)*SYZ(I,KP,L,JV)
     +             +3.*(ALF(I,K)*SZ(I,KP,L,JV)-ALF1(I,K)*FZ(I,K,JV))
         SSX (I,KP,L,JV)=SSX (I,KP,L,JV)+FX (I,K,JV)
         SZ (I,KP,L,JV)=SZ (I,KP,L,JV)+FZ (I,K,JV)
         SSXX(I,KP,L,JV)=SSXX(I,KP,L,JV)+FXX(I,K,JV)
         SSXZ(I,KP,L,JV)=SSXZ(I,KP,L,JV)+FXZ(I,K,JV)
         SZZ(I,KP,L,JV)=SZZ(I,KP,L,JV)+FZZ(I,K,JV)
C
         ENDIF
C
 330  CONTINUE
 33   CONTINUE
c     print*,'ap ADVYP 33'
C
C  traitement special pour le pole Sud (idem pole Nord)
C
      K=LAT
C
      SM0=0.
      DO 40 JV=1,NTRA
         S00(JV)=0.
 40   CONTINUE
C
      DO 41 I=1,LON
C
         IF(VGRI(I,K,L).GE.0.) THEN
           FM(I,K)=VGRI(I,K,L)*DTY
           ALF(I,K)=FM(I,K)/SM(I,K,L)
           SM(I,K,L)=SM(I,K,L)-FM(I,K)
           SM0=SM0+FM(I,K)
         ENDIF
C
         ALFQ(I,K)=ALF(I,K)*ALF(I,K)
         ALF1(I,K)=1.-ALF(I,K)
         ALF1Q(I,K)=ALF1(I,K)*ALF1(I,K)
         ALF2(I,K)=ALF1(I,K)-ALF(I,K)
         ALF3(I,K)=ALF(I,K)*ALFQ(I,K)
         ALF4(I,K)=ALF1(I,K)*ALF1Q(I,K)
C
 41   CONTINUE
c     print*,'ap ADVYP 41'
C
      DO 42 JV=1,NTRA
      DO 420 I=1,LON
C
         IF(VGRI(I,K,L).GE.0.) THEN
           F0 (I,K,JV)=ALF(I,K)* ( S0(I,K,L,JV)+ALF1(I,K)*
     +             ( SY(I,K,L,JV)+ALF2(I,K)*SYY(I,K,L,JV) ) )
           S00(JV)=S00(JV)+F0(I,K,JV)
C
           S0 (I,K,L,JV)=S0 (I,K,L,JV)-F0 (I,K,JV)
           SY (I,K,L,JV)=ALF1Q(I,K)*
     +                  (SY(I,K,L,JV)-3.*ALF(I,K)*SYY(I,K,L,JV))
           SYY(I,K,L,JV)=ALF4 (I,K)*SYY(I,K,L,JV)
      SSX (I,K,L,JV)=ALF1(I,K)*(SSX(I,K,L,JV)-ALF(I,K)*SSXY(I,K,L,JV))
      SZ (I,K,L,JV)=ALF1(I,K)*(SZ(I,K,L,JV)-ALF(I,K)*SYZ(I,K,L,JV))
           SSXX(I,K,L,JV)=ALF1 (I,K)*SSXX(I,K,L,JV)
           SSXZ(I,K,L,JV)=ALF1 (I,K)*SSXZ(I,K,L,JV)
           SZZ(I,K,L,JV)=ALF1 (I,K)*SZZ(I,K,L,JV)
           SSXY(I,K,L,JV)=ALF1Q(I,K)*SSXY(I,K,L,JV)
           SYZ(I,K,L,JV)=ALF1Q(I,K)*SYZ(I,K,L,JV)
         ENDIF
C
 420  CONTINUE
 42   CONTINUE
c     print*,'ap ADVYP 42'
C
      DO 43 I=1,LON
         IF(VGRI(I,K,L).LT.0.) THEN
           FM(I,K)=-VGRI(I,K,L)*DTY
           ALF(I,K)=FM(I,K)/SM0
         ENDIF
 43   CONTINUE
c     print*,'ap ADVYP 43'
C
      DO 44 JV=1,NTRA
      DO 440 I=1,LON
         IF(VGRI(I,K,L).LT.0.) THEN
           F0(I,K,JV)=ALF(I,K)*S00(JV)
         ENDIF
 440  CONTINUE
 44   CONTINUE
C
C  puts the temporary moments Fi into appropriate neighboring boxes
C
      DO 45 I=1,LON
C
         IF(VGRI(I,K,L).LT.0.) THEN
           SM(I,K,L)=SM(I,K,L)+FM(I,K)
           ALF(I,K)=FM(I,K)/SM(I,K,L)
         ENDIF
C
         ALFQ(I,K)=ALF(I,K)*ALF(I,K)
         ALF1(I,K)=1.-ALF(I,K)
         ALF1Q(I,K)=ALF1(I,K)*ALF1(I,K)
         ALF2(I,K)=ALF1(I,K)-ALF(I,K)
         ALF3(I,K)=ALF1(I,K)*ALF(I,K)
C
 45   CONTINUE
c     print*,'ap ADVYP 45'
C
      DO 46 JV=1,NTRA
      DO 460 I=1,LON
C
         IF(VGRI(I,K,L).LT.0.) THEN
C
         TEMPTM=-ALF(I,K)*S0(I,K,L,JV)+ALF1(I,K)*F0(I,K,JV)
         S0 (I,K,L,JV)=S0(I,K,L,JV)+F0(I,K,JV)
         SYY(I,K,L,JV)=ALF1Q(I,K)*SYY(I,K,L,JV)
     +           +5.*(-ALF3 (I,K)*SY (I,K,L,JV)+ALF2(I,K)*TEMPTM )
         SY (I,K,L,JV)=ALF1(I,K)*SY (I,K,L,JV)+3.*TEMPTM
      SSXY(I,K,L,JV)=ALF1(I,K)*SSXY(I,K,L,JV)-3.*ALF(I,K)*SSX(I,K,L,JV)
      SYZ(I,K,L,JV)=ALF1(I,K)*SYZ(I,K,L,JV)-3.*ALF(I,K)*SZ(I,K,L,JV)
C
         ENDIF
C
 460  CONTINUE
 46   CONTINUE
c     print*,'ap ADVYP 46'
C
 1    CONTINUE

c--------------------------------------------------
C     bouclage cyclique horizontal .
     
      DO l = 1,llm
         DO jv = 1,ntra
            DO j = 1,jjp1
               SM(iip1,j,l) = SM(1,j,l)
               S0(iip1,j,l,jv) = S0(1,j,l,jv)
               SSX(iip1,j,l,jv) = SSX(1,j,l,jv)   
               SY(iip1,j,l,jv) = SY(1,j,l,jv)
               SZ(iip1,j,l,jv) = SZ(1,j,l,jv)
            END DO
         END DO
      END DO

c -------------------------------------------------------------------
C *** Test  negativite:

c      DO jv = 1,ntra
c       DO l = 1,llm
c         DO j = 1,jjp1
c           DO i = 1,iip1
c              IF (s0( i,j,l,jv ).lt.0.) THEN
c                 PRINT*, '------ S0 < 0 en FIN ADVYP ---'
c                 PRINT*, 'S0(',i,j,l,jv,')=', S0(i,j,l,jv)
cc                 STOP
c              ENDIF
c           ENDDO
c         ENDDO
c       ENDDO
c      ENDDO
 
   
c -------------------------------------------------------------------
C *** Test : diag de la qtite totale de traceur dans
C            l'atmosphere avant l'advection en Y
 
       DO l = 1,llm
         DO j = 1,jjp1
           DO i = 1,iim
              sqf = sqf + S0(i,j,l,ntra)
           END DO
         END DO
       END DO
      PRINT*,'---------- DIAG DANS ADVY - SORTIE --------'
      PRINT*,'sqf=',sqf
c     print*,'ap ADVYP fin'

c-----------------------------------------------------------------
C
      RETURN
      END












