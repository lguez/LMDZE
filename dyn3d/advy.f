!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/advy.F,v 1.1.1.1 2004/05/19 12:53:06 lmdzadmin Exp $
!
      SUBROUTINE advy(limit,dty,pbarv,sm,s0,sx,sy,sz)
      use dimens_m
      use paramet_m
      use comconst
      use disvert_m
      use comgeom
      IMPLICIT NONE

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                C
C  first-order moments (SOM) advection of tracer in Y direction  C
C                                                                C
C  Source : Pascal Simon ( Meteo, CNRM )			 C
C  Adaptation : A.A. (LGGE) 					 C
C  Derniere Modif : 15/12/94 LAST
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
      PARAMETER (ntra = 1)

      REAL dty
      REAL, intent(in):: pbarv ( iip1,jjm, llm )

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

      REAL VGRI(iip1,0:jjp1,llm)

C  Rem : UGRI et WGRI ne sont pas utilises dans 
C  cette subroutine ( advection en y uniquement )
C  Rem 2 :le dimensionnement de VGRI depend de celui de pbarv
C
C  the moments F are similarly defined and used as temporary
C  storage for portions of the grid boxes in transit
C
      REAL F0(iim,0:jjp1,ntra),FM(iim,0:jjp1)
      REAL FX(iim,jjm,ntra),FY(iim,jjm,ntra)
      REAL FZ(iim,jjm,ntra)
      REAL S00(ntra)
      REAL SM0             ! Just temporal variable
C
C  work arrays
C
      REAL ALF(iim,0:jjp1),ALF1(iim,0:jjp1)
      REAL ALFQ(iim,0:jjp1),ALF1Q(iim,0:jjp1)
      REAL TEMPTM          ! Just temporal variable
c
C  Special pour poles 
c
      REAL sbms,sfms,sfzs,sbmn,sfmn,sfzn
      REAL sns0(ntra),snsz(ntra),snsm
      REAL s1v(llm),slatv(llm)
      REAL qy1(iim,llm,ntra),qylat(iim,llm,ntra)
      REAL cx1(llm,ntra), cxLAT(llm,ntra)
      REAL cy1(llm,ntra), cyLAT(llm,ntra)
      REAL z1(iim), zcos(iim), zsin(iim)
      real smpn,smps,s0pn,s0ps
      REAL SSUM
      EXTERNAL SSUM
C
      REAL sqi,sqf
      LOGICAL LIMIT

      lon = iim         ! rem : Il est possible qu'un pbl. arrive ici
      lat = jjp1        ! a cause des dim. differentes entre les
      niv=llm

C
C  the moments Fi are used as temporary storage for
C  portions of the grid boxes in transit at the current level
C
C  work arrays
C

      DO l = 1,llm
         DO j = 1,jjm
            DO i = 1,iip1  
            vgri (i,j,llm+1-l)=-1.*pbarv(i,j,l)  
            enddo
         enddo
         do i=1,iip1
             vgri(i,0,l) = 0.
             vgri(i,jjp1,l) = 0.
         enddo
      enddo

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
         sy(I,K,L,JV)=SIGN(AMIN1(AMAX1(S0(I,K,L,JV),0.),
     +                           ABS(sy(I,K,L,JV))),sy(I,K,L,JV))
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
C
 21   CONTINUE
C
      DO 22 JV=1,NTRA
      DO 220 I=1,LON
C
         IF(VGRI(I,0,L).LE.0.) THEN
C
           F0(I,0,JV)=ALF(I,0)*
     +               ( S0(I,1,L,JV)-ALF1(I,0)*sy(I,1,L,JV) )
C
           S00(JV)=S00(JV)+F0(I,0,JV)
           S0(I,1,L,JV)=S0(I,1,L,JV)-F0(I,0,JV)
           sy(I,1,L,JV)=ALF1Q(I,0)*sy(I,1,L,JV)
           sx(I,1,L,JV)=ALF1 (I,0)*sx(I,1,L,JV)
           sz(I,1,L,JV)=ALF1 (I,0)*sz(I,1,L,JV)
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
      DO 25 I=1,LON
C
         IF(VGRI(I,0,L).GT.0.) THEN
           SM(I,1,L)=SM(I,1,L)+FM(I,0)
           ALF(I,0)=FM(I,0)/SM(I,1,L)
         ENDIF
C
         ALF1(I,0)=1.-ALF(I,0)
C
 25   CONTINUE
C
      DO 26 JV=1,NTRA
      DO 260 I=1,LON
C
         IF(VGRI(I,0,L).GT.0.) THEN
C
         TEMPTM=ALF(I,0)*S0(I,1,L,JV)-ALF1(I,0)*F0(I,0,JV)
         S0(I,1,L,JV)=S0(I,1,L,JV)+F0(I,0,JV)
         sy(I,1,L,JV)=ALF1(I,0)*sy(I,1,L,JV)+3.*TEMPTM
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
C
 300  CONTINUE
 30   CONTINUE
C
      DO 31 JV=1,NTRA
      DO 31 K=1,LAT-1
      KP=K+1
      DO 310 I=1,LON
C
         IF(VGRI(I,K,L).LT.0.) THEN
C
           F0(I,K,JV)=ALF (I,K)*
     +                ( S0(I,KP,L,JV)-ALF1(I,K)*sy(I,KP,L,JV) )
           FY(I,K,JV)=ALFQ(I,K)*sy(I,KP,L,JV)
           FX(I,K,JV)=ALF (I,K)*sx(I,KP,L,JV)
           FZ(I,K,JV)=ALF (I,K)*sz(I,KP,L,JV)
C
           S0(I,KP,L,JV)=S0(I,KP,L,JV)-F0(I,K,JV)
           sy(I,KP,L,JV)=ALF1Q(I,K)*sy(I,KP,L,JV)
           sx(I,KP,L,JV)=sx(I,KP,L,JV)-FX(I,K,JV)
           sz(I,KP,L,JV)=sz(I,KP,L,JV)-FZ(I,K,JV)
C
         ELSE
C
           F0(I,K,JV)=ALF (I,K)*
     +               ( S0(I,K,L,JV)+ALF1(I,K)*sy(I,K,L,JV) )
           FY(I,K,JV)=ALFQ(I,K)*sy(I,K,L,JV)
           FX(I,K,JV)=ALF(I,K)*sx(I,K,L,JV)
           FZ(I,K,JV)=ALF(I,K)*sz(I,K,L,JV)
C
           S0(I,K,L,JV)=S0(I,K,L,JV)-F0(I,K,JV)
           sy(I,K,L,JV)=ALF1Q(I,K)*sy(I,K,L,JV)
           sx(I,K,L,JV)=sx(I,K,L,JV)-FX(I,K,JV)
           sz(I,K,L,JV)=sz(I,K,L,JV)-FZ(I,K,JV)
C
         ENDIF
C
 310  CONTINUE
 31   CONTINUE
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
         ALF1(I,K)=1.-ALF(I,K)
C
 320  CONTINUE
 32   CONTINUE
C
      DO 33 JV=1,NTRA
      DO 33 K=1,LAT-1
      KP=K+1
      DO 330 I=1,LON
C
         IF(VGRI(I,K,L).LT.0.) THEN
C
         TEMPTM=-ALF(I,K)*S0(I,K,L,JV)+ALF1(I,K)*F0(I,K,JV)
         S0(I,K,L,JV)=S0(I,K,L,JV)+F0(I,K,JV)
         sy(I,K,L,JV)=ALF(I,K)*FY(I,K,JV)+ALF1(I,K)*sy(I,K,L,JV)
     +               +3.*TEMPTM
         sx(I,K,L,JV)=sx(I,K,L,JV)+FX(I,K,JV)
         sz(I,K,L,JV)=sz(I,K,L,JV)+FZ(I,K,JV)
C
         ELSE
C
         TEMPTM=ALF(I,K)*S0(I,KP,L,JV)-ALF1(I,K)*F0(I,K,JV)
         S0(I,KP,L,JV)=S0(I,KP,L,JV)+F0(I,K,JV)
         sy(I,KP,L,JV)=ALF(I,K)*FY(I,K,JV)+ALF1(I,K)*sy(I,KP,L,JV)
     +                +3.*TEMPTM
         sx(I,KP,L,JV)=sx(I,KP,L,JV)+FX(I,K,JV)
         sz(I,KP,L,JV)=sz(I,KP,L,JV)+FZ(I,K,JV)
C
         ENDIF
C
 330  CONTINUE
 33   CONTINUE
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
C
 41   CONTINUE
C
      DO 42 JV=1,NTRA
      DO 420 I=1,LON
C
         IF(VGRI(I,K,L).GE.0.) THEN
           F0 (I,K,JV)=ALF(I,K)*
     +                ( S0(I,K,L,JV)+ALF1(I,K)*sy(I,K,L,JV) )
           S00(JV)=S00(JV)+F0(I,K,JV)
C
           S0(I,K,L,JV)=S0 (I,K,L,JV)-F0 (I,K,JV)
           sy(I,K,L,JV)=ALF1Q(I,K)*sy(I,K,L,JV)
           sx(I,K,L,JV)=ALF1(I,K)*sx(I,K,L,JV)
           sz(I,K,L,JV)=ALF1(I,K)*sz(I,K,L,JV)
         ENDIF
C
 420  CONTINUE
 42   CONTINUE
C
      DO 43 I=1,LON
         IF(VGRI(I,K,L).LT.0.) THEN
           FM(I,K)=-VGRI(I,K,L)*DTY
           ALF(I,K)=FM(I,K)/SM0
         ENDIF
 43   CONTINUE
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
         ALF1(I,K)=1.-ALF(I,K)
C
 45   CONTINUE
C
      DO 46 JV=1,NTRA
      DO 460 I=1,LON
C
         IF(VGRI(I,K,L).LT.0.) THEN
C
         TEMPTM=-ALF(I,K)*S0(I,K,L,JV)+ALF1(I,K)*F0(I,K,JV)
         S0(I,K,L,JV)=S0(I,K,L,JV)+F0(I,K,JV)
         sy(I,K,L,JV)=ALF1(I,K)*sy(I,K,L,JV)+3.*TEMPTM
C
         ENDIF
C
 460  CONTINUE
 46   CONTINUE
C
 1    CONTINUE
C
      RETURN
      END

