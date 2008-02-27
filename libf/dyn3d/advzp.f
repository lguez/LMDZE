!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/advzp.F,v 1.1.1.1 2004/05/19 12:53:06 lmdzadmin Exp $
!
      SUBROUTINE ADVZP(LIMIT,DTZ,W,SM,S0,SSX,SY,SZ
     .                 ,SSXX,SSXY,SSXZ,SYY,SYZ,SZZ,ntra )

      use dimens_m
      use paramet_m
      use comconst
      use comvert
      use comgeom
      IMPLICIT NONE

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                 C
C  second-order moments (SOM) advection of tracer in Z direction  C
C                                                                 C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                 C
C  Source : Pascal Simon ( Meteo, CNRM )                          C
C  Adaptation : A.A. (LGGE)                                       C
C  Derniere Modif : 19/11/95 LAST                                 C
C                                                                 C
C  sont les arguments d'entree pour le s-pg                       C
C                                                                 C
C  argument de sortie du s-pg                                     C
C                                                                 C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Rem : Probleme aux poles il faut reecrire ce cas specifique
C        Attention au sens de l'indexation
C

C
C  parametres principaux du modele
C
C
C  Arguments :
C  ----------
C  dty : frequence fictive d'appel du transport
C  parbu,pbarv : flux de masse en x et y en Pa.m2.s-1
c
        INTEGER lon,lat,niv
        INTEGER i,j,jv,k,kp,l,lp
        INTEGER ntra
c        PARAMETER (ntra = 1)
c
        REAL dtz
        REAL w ( iip1,jjp1,llm )
c
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
C
C  mass fluxes across the boundaries (UGRI,VGRI,WGRI)
C  mass fluxes in kg
C  declaration :
C
      REAL WGRI(iip1,jjp1,0:llm)

C Rem : UGRI et VGRI ne sont pas utilises dans
C  cette subroutine ( advection en z uniquement )
C  Rem 2 :le dimensionnement de VGRI depend de celui de pbarv
C         attention a celui de WGRI
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
      REAL F0(iim,llm,ntra),FM(iim,llm)
      REAL FX(iim,llm,ntra),FY(iim,llm,ntra)
      REAL FZ(iim,llm,ntra)
      REAL FXX(iim,llm,ntra),FXY(iim,llm,ntra)
      REAL FXZ(iim,llm,ntra),FYY(iim,llm,ntra)
      REAL FYZ(iim,llm,ntra),FZZ(iim,llm,ntra)
      REAL S00(ntra)
      REAL SM0             ! Just temporal variable
C
C  work arrays
C
      REAL ALF(iim),ALF1(iim)
      REAL ALFQ(iim),ALF1Q(iim)
      REAL ALF2(iim),ALF3(iim)
      REAL ALF4(iim)
      REAL TEMPTM          ! Just temporal variable
      REAL SLPMAX,S1MAX,S1NEW,S2NEW
c
      REAL sqi,sqf
      LOGICAL LIMIT

      lon = iim         ! rem : Il est possible qu'un pbl. arrive ici
      lat = jjp1        ! a cause des dim. differentes entre les
      niv = llm         !       tab. S et VGRI 
                    
c-----------------------------------------------------------------
C *** Test : diag de la qtite totale de traceur dans
C            l'atmosphere avant l'advection en Y
c 
      sqi = 0.
      sqf = 0.
c
      DO l = 1,llm
         DO j = 1,jjp1
           DO i = 1,iim
              sqi = sqi + S0(i,j,l,ntra)
           END DO
         END DO
      END DO
      PRINT*,'---------- DIAG DANS ADVZP - ENTREE --------'
      PRINT*,'sqi=',sqi

c-----------------------------------------------------------------
C  Interface : adaptation nouveau modele
C  -------------------------------------
C
C  Conversion des flux de masses en kg

      DO 500 l = 1,llm
         DO 500 j = 1,jjp1
            DO 500 i = 1,iip1  
            wgri (i,j,llm+1-l) = w (i,j,l)  
  500 CONTINUE
      do j=1,jjp1
         do i=1,iip1
            wgri(i,j,0)=0.
         enddo
      enddo
c
cAA rem : Je ne suis pas sur du signe  
cAA       Je ne suis pas sur pour le 0:llm
c
c-----------------------------------------------------------------
C---------------------- START HERE -------------------------------
C
C  boucle sur les latitudes
C
      DO 1 K=1,LAT
C
C  place limits on appropriate moments before transport
C      (if flux-limiting is to be applied)
C
      IF(.NOT.LIMIT) GO TO 101
C
      DO 10 JV=1,NTRA
      DO 10 L=1,NIV
         DO 100 I=1,LON
            IF(S0(I,K,L,JV).GT.0.) THEN
              SLPMAX=S0(I,K,L,JV)
              S1MAX =1.5*SLPMAX
              S1NEW =AMIN1(S1MAX,AMAX1(-S1MAX,SZ(I,K,L,JV)))
              S2NEW =AMIN1( 2.*SLPMAX-ABS(S1NEW)/3. ,
     +                     AMAX1(ABS(S1NEW)-SLPMAX,SZZ(I,K,L,JV)) )
              SZ (I,K,L,JV)=S1NEW
              SZZ(I,K,L,JV)=S2NEW
              SSXZ(I,K,L,JV)=AMIN1(SLPMAX,AMAX1(-SLPMAX,SSXZ(I,K,L,JV)))
              SYZ(I,K,L,JV)=AMIN1(SLPMAX,AMAX1(-SLPMAX,SYZ(I,K,L,JV)))
            ELSE
              SZ (I,K,L,JV)=0.
              SZZ(I,K,L,JV)=0.
              SSXZ(I,K,L,JV)=0.
              SYZ(I,K,L,JV)=0.
            ENDIF
 100     CONTINUE
 10   CONTINUE
C
 101  CONTINUE
C
C  boucle sur les niveaux intercouches de 1 a NIV-1
C   (flux nul au sommet L=0 et a la base L=NIV)
C
C  calculate flux and moments between adjacent boxes
C     (flux from LP to L if WGRI(L).lt.0, from L to LP if WGRI(L).gt.0)
C  1- create temporary moments/masses for partial boxes in transit
C  2- reajusts moments remaining in the box
C
      DO 11 L=1,NIV-1
      LP=L+1
C
      DO 110 I=1,LON
C
         IF(WGRI(I,K,L).LT.0.) THEN
           FM(I,L)=-WGRI(I,K,L)*DTZ
           ALF(I)=FM(I,L)/SM(I,K,LP)
           SM(I,K,LP)=SM(I,K,LP)-FM(I,L)
         ELSE
           FM(I,L)=WGRI(I,K,L)*DTZ
           ALF(I)=FM(I,L)/SM(I,K,L)
           SM(I,K,L)=SM(I,K,L)-FM(I,L)
         ENDIF
C
         ALFQ (I)=ALF(I)*ALF(I)
         ALF1 (I)=1.-ALF(I)
         ALF1Q(I)=ALF1(I)*ALF1(I)
         ALF2 (I)=ALF1(I)-ALF(I)
         ALF3 (I)=ALF(I)*ALFQ(I)
         ALF4 (I)=ALF1(I)*ALF1Q(I)
C
 110  CONTINUE
C
      DO 111 JV=1,NTRA
      DO 1110 I=1,LON
C
         IF(WGRI(I,K,L).LT.0.) THEN
C
           F0 (I,L,JV)=ALF (I)* ( S0(I,K,LP,JV)-ALF1(I)*
     +          ( SZ(I,K,LP,JV)-ALF2(I)*SZZ(I,K,LP,JV) ) )
           FZ (I,L,JV)=ALFQ(I)*(SZ(I,K,LP,JV)-3.*ALF1(I)*SZZ(I,K,LP,JV))
           FZZ(I,L,JV)=ALF3(I)*SZZ(I,K,LP,JV)
           FXZ(I,L,JV)=ALFQ(I)*SSXZ(I,K,LP,JV)
           FYZ(I,L,JV)=ALFQ(I)*SYZ(I,K,LP,JV)
           FX (I,L,JV)=ALF (I)*(SSX(I,K,LP,JV)-ALF1(I)*SSXZ(I,K,LP,JV))
           FY (I,L,JV)=ALF (I)*(SY(I,K,LP,JV)-ALF1(I)*SYZ(I,K,LP,JV))
           FXX(I,L,JV)=ALF (I)*SSXX(I,K,LP,JV)
           FXY(I,L,JV)=ALF (I)*SSXY(I,K,LP,JV)
           FYY(I,L,JV)=ALF (I)*SYY(I,K,LP,JV)
C
           S0 (I,K,LP,JV)=S0 (I,K,LP,JV)-F0 (I,L,JV)
           SZ (I,K,LP,JV)=ALF1Q(I)
     +                   *(SZ(I,K,LP,JV)+3.*ALF(I)*SZZ(I,K,LP,JV))
           SZZ(I,K,LP,JV)=ALF4 (I)*SZZ(I,K,LP,JV)
           SSXZ(I,K,LP,JV)=ALF1Q(I)*SSXZ(I,K,LP,JV)
           SYZ(I,K,LP,JV)=ALF1Q(I)*SYZ(I,K,LP,JV)
           SSX (I,K,LP,JV)=SSX (I,K,LP,JV)-FX (I,L,JV)
           SY (I,K,LP,JV)=SY (I,K,LP,JV)-FY (I,L,JV)
           SSXX(I,K,LP,JV)=SSXX(I,K,LP,JV)-FXX(I,L,JV)
           SSXY(I,K,LP,JV)=SSXY(I,K,LP,JV)-FXY(I,L,JV)
           SYY(I,K,LP,JV)=SYY(I,K,LP,JV)-FYY(I,L,JV)
C
         ELSE
C
           F0 (I,L,JV)=ALF (I)*(S0(I,K,L,JV)
     +           +ALF1(I) * (SZ(I,K,L,JV)+ALF2(I)*SZZ(I,K,L,JV)) )
           FZ (I,L,JV)=ALFQ(I)*(SZ(I,K,L,JV)+3.*ALF1(I)*SZZ(I,K,L,JV))
           FZZ(I,L,JV)=ALF3(I)*SZZ(I,K,L,JV)
           FXZ(I,L,JV)=ALFQ(I)*SSXZ(I,K,L,JV)
           FYZ(I,L,JV)=ALFQ(I)*SYZ(I,K,L,JV)
           FX (I,L,JV)=ALF (I)*(SSX(I,K,L,JV)+ALF1(I)*SSXZ(I,K,L,JV))
           FY (I,L,JV)=ALF (I)*(SY(I,K,L,JV)+ALF1(I)*SYZ(I,K,L,JV))
           FXX(I,L,JV)=ALF (I)*SSXX(I,K,L,JV)
           FXY(I,L,JV)=ALF (I)*SSXY(I,K,L,JV)
           FYY(I,L,JV)=ALF (I)*SYY(I,K,L,JV)
C
           S0 (I,K,L,JV)=S0 (I,K,L,JV)-F0(I,L,JV)
           SZ (I,K,L,JV)=ALF1Q(I)*(SZ(I,K,L,JV)-3.*ALF(I)*SZZ(I,K,L,JV))
           SZZ(I,K,L,JV)=ALF4 (I)*SZZ(I,K,L,JV)
           SSXZ(I,K,L,JV)=ALF1Q(I)*SSXZ(I,K,L,JV)
           SYZ(I,K,L,JV)=ALF1Q(I)*SYZ(I,K,L,JV)
           SSX (I,K,L,JV)=SSX (I,K,L,JV)-FX (I,L,JV)
           SY (I,K,L,JV)=SY (I,K,L,JV)-FY (I,L,JV)
           SSXX(I,K,L,JV)=SSXX(I,K,L,JV)-FXX(I,L,JV)
           SSXY(I,K,L,JV)=SSXY(I,K,L,JV)-FXY(I,L,JV)
           SYY(I,K,L,JV)=SYY(I,K,L,JV)-FYY(I,L,JV)
C
         ENDIF
C
 1110 CONTINUE
 111  CONTINUE
C
 11   CONTINUE
C
C  puts the temporary moments Fi into appropriate neighboring boxes
C
      DO 12 L=1,NIV-1
      LP=L+1
C
      DO 120 I=1,LON
C
         IF(WGRI(I,K,L).LT.0.) THEN
           SM(I,K,L)=SM(I,K,L)+FM(I,L)
           ALF(I)=FM(I,L)/SM(I,K,L)
         ELSE
           SM(I,K,LP)=SM(I,K,LP)+FM(I,L)
           ALF(I)=FM(I,L)/SM(I,K,LP)
         ENDIF
C
         ALF1(I)=1.-ALF(I)
         ALFQ(I)=ALF(I)*ALF(I)
         ALF1Q(I)=ALF1(I)*ALF1(I)
         ALF2(I)=ALF(I)*ALF1(I)
         ALF3(I)=ALF1(I)-ALF(I)
C
 120  CONTINUE
C
      DO 121 JV=1,NTRA
      DO 1210 I=1,LON
C
         IF(WGRI(I,K,L).LT.0.) THEN
C
           TEMPTM=-ALF(I)*S0(I,K,L,JV)+ALF1(I)*F0(I,L,JV)
           S0 (I,K,L,JV)=S0(I,K,L,JV)+F0(I,L,JV)
           SZZ(I,K,L,JV)=ALFQ(I)*FZZ(I,L,JV)+ALF1Q(I)*SZZ(I,K,L,JV)
     +        +5.*( ALF2(I)*(FZ(I,L,JV)-SZ(I,K,L,JV))+ALF3(I)*TEMPTM )
           SZ (I,K,L,JV)=ALF (I)*FZ (I,L,JV)+ALF1 (I)*SZ (I,K,L,JV)
     +                  +3.*TEMPTM
           SSXZ(I,K,L,JV)=ALF (I)*FXZ(I,L,JV)+ALF1 (I)*SSXZ(I,K,L,JV)
     +              +3.*(ALF1(I)*FX (I,L,JV)-ALF  (I)*SSX (I,K,L,JV))
           SYZ(I,K,L,JV)=ALF (I)*FYZ(I,L,JV)+ALF1 (I)*SYZ(I,K,L,JV)
     +              +3.*(ALF1(I)*FY (I,L,JV)-ALF  (I)*SY (I,K,L,JV))
           SSX (I,K,L,JV)=SSX (I,K,L,JV)+FX (I,L,JV)
           SY (I,K,L,JV)=SY (I,K,L,JV)+FY (I,L,JV)
           SSXX(I,K,L,JV)=SSXX(I,K,L,JV)+FXX(I,L,JV)
           SSXY(I,K,L,JV)=SSXY(I,K,L,JV)+FXY(I,L,JV)
           SYY(I,K,L,JV)=SYY(I,K,L,JV)+FYY(I,L,JV)
C
         ELSE
C
           TEMPTM=ALF(I)*S0(I,K,LP,JV)-ALF1(I)*F0(I,L,JV)
           S0 (I,K,LP,JV)=S0(I,K,LP,JV)+F0(I,L,JV)
           SZZ(I,K,LP,JV)=ALFQ(I)*FZZ(I,L,JV)+ALF1Q(I)*SZZ(I,K,LP,JV)
     +        +5.*( ALF2(I)*(SZ(I,K,LP,JV)-FZ(I,L,JV))-ALF3(I)*TEMPTM )
           SZ (I,K,LP,JV)=ALF (I)*FZ(I,L,JV)+ALF1(I)*SZ(I,K,LP,JV)
     +                   +3.*TEMPTM
           SSXZ(I,K,LP,JV)=ALF(I)*FXZ(I,L,JV)+ALF1(I)*SSXZ(I,K,LP,JV)
     +                   +3.*(ALF(I)*SSX(I,K,LP,JV)-ALF1(I)*FX(I,L,JV))
           SYZ(I,K,LP,JV)=ALF(I)*FYZ(I,L,JV)+ALF1(I)*SYZ(I,K,LP,JV)
     +                   +3.*(ALF(I)*SY(I,K,LP,JV)-ALF1(I)*FY(I,L,JV))
           SSX (I,K,LP,JV)=SSX (I,K,LP,JV)+FX (I,L,JV)
           SY (I,K,LP,JV)=SY (I,K,LP,JV)+FY (I,L,JV)
           SSXX(I,K,LP,JV)=SSXX(I,K,LP,JV)+FXX(I,L,JV)
           SSXY(I,K,LP,JV)=SSXY(I,K,LP,JV)+FXY(I,L,JV)
           SYY(I,K,LP,JV)=SYY(I,K,LP,JV)+FYY(I,L,JV)
C
         ENDIF
C
 1210 CONTINUE
 121  CONTINUE
C
 12   CONTINUE
C
C  fin de la boucle principale sur les latitudes
C
 1    CONTINUE
C
      DO l = 1,llm
      DO j = 1,jjp1
          SM(iip1,j,l) = SM(1,j,l)
	  S0(iip1,j,l,ntra) = S0(1,j,l,ntra)
          SSX(iip1,j,l,ntra) = SSX(1,j,l,ntra)
	  SY(iip1,j,l,ntra) = SY(1,j,l,ntra)
          SZ(iip1,j,l,ntra) = SZ(1,j,l,ntra)
      ENDDO
      ENDDO
c										C-------------------------------------------------------------
C *** Test : diag de la qqtite totale de tarceur
C            dans l'atmosphere avant l'advection en z
       DO l = 1,llm
       DO j = 1,jjp1
       DO i = 1,iim
          sqf = sqf + S0(i,j,l,ntra)
       ENDDO
       ENDDO
       ENDDO
       PRINT*,'-------- DIAG DANS ADVZ - SORTIE ---------'
       PRINT*,'sqf=', sqf

      RETURN
      END
