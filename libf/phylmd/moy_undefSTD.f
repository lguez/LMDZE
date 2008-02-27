      SUBROUTINE moy_undefSTD(nlevSTD,itap,
     $           ecrit_day,ecrit_mth,ecrit_hf2mth,
     $           tnondef,tsumSTD)
      use dimens_m
      use dimphy
      IMPLICIT none
c moyenne d'1 var a 1 niveau de pression
c====================================================================
c
c I. Musat : 09.2004
c
c Moyenne - a des frequences differentes - des valeurs bien definies
c         (.NE.1.E+20) des variables interpolees a un niveau de
c         pression.
c 1) les variables de type "day" (nout=1) ou "mth" (nout=2) sont sommees
c    tous les pas de temps de la physique
c
c 2) les variables de type "NMC" (nout=3) sont calculees a partir
c    des valeurs instantannees toutes les 6 heures
c
c
c NB: mettre "inst(X)" dans le write_histXXX.h !
c====================================================================
      integer jjmp1
      parameter (jjmp1=jjm+1-1/jjm)
c
c variables Input
      INTEGER nlevSTD, klevSTD, itap
      PARAMETER(klevSTD=17)
      INTEGER ecrit_day, ecrit_mth, ecrit_hf2mth
c
c variables locales
      INTEGER i, k, nout
      PARAMETER(nout=3) !nout=1 day/nout=2 mth/nout=3 NMC
c
c variables Output
      REAL tnondef(klon,klevSTD,nout)
      REAL tsumSTD(klon,klevSTD,nout)
c
c calcul 1 fois par jour
c
      IF(MOD(itap,ecrit_day).EQ.0) THEN
       DO k=1, nlevSTD
        DO i=1, klon
         IF(tnondef(i,k,1).NE.FLOAT(ecrit_day)) THEN
          tsumSTD(i,k,1)=tsumSTD(i,k,1)/
     $    (FLOAT(ecrit_day)-tnondef(i,k,1))
         ELSE
          tsumSTD(i,k,1)=1.E+20
         ENDIF !tnondef
        ENDDO !i
       ENDDO !k
      ENDIF !MOD(itap,ecrit_day).EQ.0
c
c calcul 1 fois par mois
c
      IF(MOD(itap,ecrit_mth).EQ.0) THEN
       DO k=1, nlevSTD
        DO i=1, klon
         IF(tnondef(i,k,2).NE.FLOAT(ecrit_mth)) THEN
          tsumSTD(i,k,2)=tsumSTD(i,k,2)/
     $    (FLOAT(ecrit_mth)-tnondef(i,k,2))
         ELSE
          tsumSTD(i,k,2)=1.E+20
         ENDIF !tnondef
c
         IF(tnondef(i,k,3).NE.FLOAT(ecrit_hf2mth)) THEN
          tsumSTD(i,k,3)=tsumSTD(i,k,3)/
     $    (FLOAT(ecrit_hf2mth)-tnondef(i,k,3))
         ELSE
          tsumSTD(i,k,3)=1.E+20
         ENDIF !tnondef
c
        ENDDO !i
       ENDDO !k
      ENDIF !MOD
c
      RETURN
      END  
