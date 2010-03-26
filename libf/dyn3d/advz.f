!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/advz.F,v 1.2 2005/05/25 13:10:09 fairhead Exp $
!
      SUBROUTINE advz(limit,dtz,w,sm,s0,sx,sy,sz)
      use dimens_m
      use paramet_m
      use comconst
      use comvert
      IMPLICIT NONE

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                C
C  first-order moments (FOM) advection of tracer in Z direction  C
C                                                                C
C  Source : Pascal Simon (Meteo,CNRM)                            C
C  Adaptation : A.Armengaud (LGGE) juin 94                       C
C                                                                C
C                                                                C
C  sont des arguments d'entree pour le s-pg...                   C
C                                                                C
C  dq est l'argument de sortie pour le s-pg                      C
C								 C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  parametres principaux du modele
C

C  Arguments :
C  -----------
C  dtz : frequence fictive d'appel du transport 
C  w : flux de masse en z en Pa.m2.s-1

      INTEGER ntra
      PARAMETER (ntra = 1)

      REAL, intent(in):: dtz
      REAL w ( iip1,jjp1,llm )
    
C  moments: SM  total mass in each grid box
C           S0  mass of tracer in each grid box
C           Si  1rst order moment in i direction
C
      REAL SM(iip1,jjp1,llm)
     +    ,S0(iip1,jjp1,llm,ntra)
      REAL sx(iip1,jjp1,llm,ntra)
     +    ,sy(iip1,jjp1,llm,ntra)
     +    ,sz(iip1,jjp1,llm,ntra)


C  Local :
C  ------- 

C  mass fluxes across the boundaries (UGRI,VGRI,WGRI)
C  mass fluxes in kg
C  declaration :

      REAL WGRI(iip1,jjp1,0:llm)

C
C  the moments F are used as temporary  storage for 
C  portions of grid boxes in transit at the current latitude
C
      REAL FM(iim,llm)
      REAL F0(iim,llm,ntra),FX(iim,llm,ntra)
      REAL FY(iim,llm,ntra),FZ(iim,llm,ntra)
C
C  work arrays
C
      REAL ALF(iim),ALF1(iim),ALFQ(iim),ALF1Q(iim)
      REAL TEMPTM            ! Just temporal variable
      REAL sqi,sqf
C
      LOGICAL LIMIT
      INTEGER lon,lat,niv
      INTEGER i,j,jv,k,l,lp

      lon = iim
      lat = jjp1
      niv = llm 

C *** Test de passage d'arguments ******
 
c     DO 399 l = 1, llm
c     DO 399 j = 1, jjp1
c     DO 399 i = 1, iip1
c        IF (S0(i,j,l,ntra) .lt. 0. ) THEN
c           PRINT*,'S0(',i,j,l,')=',S0(i,j,l,ntra)
c           print*, 'sx(',i,j,l,')=',sx(i,j,l,ntra)
c           print*, 'sy(',i,j,l,')=',sy(i,j,l,ntra)
c           print*, 'sz(',i,j,l,')=',sz(i,j,l,ntra)
c           PRINT*, 'AIE !! debut ADVZ - pbl arg. passage dans ADVZ'
c            STOP
c        ENDIF
  399 CONTINUE

C-----------------------------------------------------------------
C *** Test : diag de la qqtite totale de traceur 
C            dans l'atmosphere avant l'advection en z
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
      PRINT*,'-------- DIAG DANS ADVZ - ENTREE ---------'
      PRINT*,'sqi=',sqi

C-----------------------------------------------------------------
C  Interface : adaptation nouveau modele
C  -------------------------------------
C
C  Conversion du flux de masse en kg.s-1

      DO 500 l = 1,llm
         DO 500 j = 1,jjp1
            DO 500 i = 1,iip1  
c            wgri (i,j,llm+1-l) =  w (i,j,l) / g 
               wgri (i,j,llm+1-l) =  w (i,j,l) 
c             wgri (i,j,0) = 0.                ! a detruire ult.
c             wgri (i,j,l) = 0.1               !    w (i,j,l) 
c             wgri (i,j,llm) = 0.              ! a detruire ult.
  500 CONTINUE
         DO  j = 1,jjp1
            DO i = 1,iip1  
               wgri(i,j,0)=0.
            enddo
         enddo

C-----------------------------------------------------------------
  
C  start here          
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
            sz(I,K,L,JV)=SIGN(AMIN1(AMAX1(S0(I,K,L,JV),0.),
     +                              ABS(sz(I,K,L,JV))),sz(I,K,L,JV))
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
C
 110  CONTINUE
C
      DO 111 JV=1,NTRA
      DO 1110 I=1,LON
C
         IF(WGRI(I,K,L).LT.0.) THEN
C
           F0(I,L,JV)=ALF (I)*( S0(I,K,LP,JV)-ALF1(I)*sz(I,K,LP,JV) )
           FZ(I,L,JV)=ALFQ(I)*sz(I,K,LP,JV)
           FX(I,L,JV)=ALF (I)*sx(I,K,LP,JV)
           FY(I,L,JV)=ALF (I)*sy(I,K,LP,JV)
C
           S0(I,K,LP,JV)=S0(I,K,LP,JV)-F0(I,L,JV)
           sz(I,K,LP,JV)=ALF1Q(I)*sz(I,K,LP,JV)
           sx(I,K,LP,JV)=sx(I,K,LP,JV)-FX(I,L,JV)
           sy(I,K,LP,JV)=sy(I,K,LP,JV)-FY(I,L,JV)
C
         ELSE
C
           F0(I,L,JV)=ALF (I)*(S0(I,K,L,JV)+ALF1(I)*sz(I,K,L,JV) )
           FZ(I,L,JV)=ALFQ(I)*sz(I,K,L,JV)
           FX(I,L,JV)=ALF (I)*sx(I,K,L,JV)
           FY(I,L,JV)=ALF (I)*sy(I,K,L,JV)
C
           S0(I,K,L,JV)=S0(I,K,L,JV)-F0(I,L,JV)
           sz(I,K,L,JV)=ALF1Q(I)*sz(I,K,L,JV)
           sx(I,K,L,JV)=sx(I,K,L,JV)-FX(I,L,JV)
           sy(I,K,L,JV)=sy(I,K,L,JV)-FY(I,L,JV)
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
C
 120  CONTINUE
C
      DO 121 JV=1,NTRA
      DO 1210 I=1,LON
C
         IF(WGRI(I,K,L).LT.0.) THEN
C
           TEMPTM=-ALF(I)*S0(I,K,L,JV)+ALF1(I)*F0(I,L,JV)
           S0(I,K,L,JV)=S0(I,K,L,JV)+F0(I,L,JV)
           sz(I,K,L,JV)=ALF(I)*FZ(I,L,JV)+ALF1(I)*sz(I,K,L,JV)+3.*TEMPTM
           sx(I,K,L,JV)=sx(I,K,L,JV)+FX(I,L,JV)
           sy(I,K,L,JV)=sy(I,K,L,JV)+FY(I,L,JV)
C
         ELSE
C
           TEMPTM=ALF(I)*S0(I,K,LP,JV)-ALF1(I)*F0(I,L,JV)
           S0(I,K,LP,JV)=S0(I,K,LP,JV)+F0(I,L,JV)
           sz(I,K,LP,JV)=ALF(I)*FZ(I,L,JV)+ALF1(I)*sz(I,K,LP,JV)
     +                  +3.*TEMPTM
           sx(I,K,LP,JV)=sx(I,K,LP,JV)+FX(I,L,JV)
           sy(I,K,LP,JV)=sy(I,K,LP,JV)+FY(I,L,JV)
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
C-------------------------------------------------------------
C
C ----------- AA Test en fin de ADVX ------ Controle des S*

c     DO 9999 l = 1, llm
c     DO 9999 j = 1, jjp1
c     DO 9999 i = 1, iip1
c        IF (S0(i,j,l,ntra).lt.0..and.LIMIT) THEN 
c           PRINT*, '-------------------'
c           PRINT*, 'En fin de ADVZ'
c           PRINT*,'S0(',i,j,l,')=',S0(i,j,l,ntra)
c           print*, 'sx(',i,j,l,')=',sx(i,j,l,ntra)
c           print*, 'sy(',i,j,l,')=',sy(i,j,l,ntra)
c           print*, 'sz(',i,j,l,')=',sz(i,j,l,ntra)
c           WRITE (*,*) 'On arrete !! - pbl en fin de ADVZ1'
c            STOP
c        ENDIF
 9999 CONTINUE

C *** ------------------- bouclage cyclique  en X ------------
      
c      DO l = 1,llm
c         DO j = 1,jjp1
c            SM(iip1,j,l) = SM(1,j,l)
c            S0(iip1,j,l,ntra) = S0(1,j,l,ntra)
C            sx(iip1,j,l,ntra) = sx(1,j,l,ntra)
c            sy(iip1,j,l,ntra) = sy(1,j,l,ntra)
c            sz(iip1,j,l,ntra) = sz(1,j,l,ntra)
c         ENDDO
c      ENDDO
           
C-------------------------------------------------------------
C *** Test : diag de la qqtite totale de traceur 
C            dans l'atmosphere avant l'advection en z
      DO l = 1,llm
         DO j = 1,jjp1
            DO i = 1,iim
cIM 240305            sqf = sqf + S0(i,j,l,9)
               sqf = sqf + S0(i,j,l,ntra)
            ENDDO
         ENDDO
      ENDDO
      PRINT*,'-------- DIAG DANS ADVZ - SORTIE ---------'
      PRINT*,'sqf=', sqf

C-------------------------------------------------------------
      RETURN
      END
C_______________________________________________________________
C_______________________________________________________________
