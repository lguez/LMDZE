      SUBROUTINE undefSTD(nlevSTD,itap,tlevSTD,
     $           ecrit_hf,
     $           oknondef,tnondef,tsumSTD)
      use dimens_m
      use dimphy
      IMPLICIT none
c
c====================================================================
c
c I. Musat : 09.2004
c
c Calcul * du nombre de pas de temps (FLOAT(ecrit_XXX)-tnondef)) 
c          ou la variable tlevSTD est bien definie (.NE.1.E+20), 
c et 
c        * de la somme de tlevSTD => tsumSTD
c
c nout=1 !var. journaliere "day" moyenne sur tous les pas de temps
c        ! de la physique
c nout=2 !var. mensuelle "mth" moyennee sur tous les pas de temps
c        ! de la physique
c nout=3 !var. mensuelle "NMC" moyennee toutes les 6heures
c
c
c NB: mettre "inst(X)" dans le write_histXXX.h !
c====================================================================
c
      integer jjmp1
      parameter (jjmp1=jjm+1-1/jjm)
c variables Input
      INTEGER nlevSTD, klevSTD, itap
      PARAMETER(klevSTD=17)
      INTEGER, intent(in):: ecrit_hf
c
c variables locales
      INTEGER i, k, nout
      PARAMETER(nout=3) !nout=1 : day; =2 : mth; =3 : NMC
c
c variables Output
      REAL tlevSTD(klon,klevSTD), tsumSTD(klon,klevSTD,nout)
      LOGICAL oknondef(klon,klevSTD,nout)
      REAL tnondef(klon,klevSTD,nout)
c
c calcul variables tous les pas de temps de la physique 
c
      DO k=1, nlevSTD
       DO i=1, klon
        IF(tlevSTD(i,k).EQ.1E+20) THEN
         IF(oknondef(i,k,1)) THEN          
          tnondef(i,k,1)=tnondef(i,k,1)+1.
         ENDIF !oknondef(i,k)
c
         IF(oknondef(i,k,2)) THEN          
          tnondef(i,k,2)=tnondef(i,k,2)+1.
         ENDIF !oknondef(i,k)
c
        ELSE IF(tlevSTD(i,k).NE.1E+20) THEN
         tsumSTD(i,k,1)=tsumSTD(i,k,1)+tlevSTD(i,k)
         tsumSTD(i,k,2)=tsumSTD(i,k,2)+tlevSTD(i,k)
        ENDIF 
       ENDDO !i
      ENDDO !k
c
c calcul variables toutes les 6h
c
      IF(MOD(itap,ecrit_hf).EQ.0) THEN
c
       DO k=1, nlevSTD
        DO i=1, klon
         IF(tlevSTD(i,k).EQ.1E+20) THEN
          IF(oknondef(i,k,3)) THEN          
           tnondef(i,k,3)=tnondef(i,k,3)+1.
          ENDIF !oknondef(i,k)
c
         ELSE IF(tlevSTD(i,k).NE.1E+20) THEN
         tsumSTD(i,k,3)=tsumSTD(i,k,3)+tlevSTD(i,k)
         ENDIF 
        ENDDO !i
       ENDDO !k

      ENDIF !MOD(itap,ecrit_hf).EQ.0
c
      RETURN
      END  
