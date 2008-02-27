      SUBROUTINE ini_undefSTD(nlevSTD,itap,
     $           ecrit_day,ecrit_mth,
     $           tnondef,tsumSTD)
      use dimens_m
      use dimphy
      IMPLICIT none
c
c====================================================================
c
c I. Musat : 09.2004
c
c Initialisation - a des frequences differentes : 
c
c 1) des variables moyennees sur la journee "day" ou sur le mois "mth"
c    calculees a partir des valeurs "instantannees" de la physique
c
c 2) des variables moyennes mensuelles "NMC" calculees a partir des val.
c    toutes les 6 heures
c
c nout=1 !var. journaliere "day" moyenne sur tous les pas de temps
c              ! de la physique
c nout=2 !var. mensuelle "mth" moyennee sur tous les pas de temps
c              ! de la physique
c nout=3 !var. mensuelle "NMC" moyennee toutes les 6heures
c
c
c NB: mettre "inst(X)" dans le write_histXXX.h !
c====================================================================
c
      integer jjmp1
      parameter (jjmp1=jjm+1-1/jjm)
c variables Input/Output
      INTEGER, intent(in):: nlevSTD
      integer klevSTD, itap
      PARAMETER(klevSTD=17)
      INTEGER ecrit_day,ecrit_mth
c
c variables locales
      INTEGER i, k, nout
      PARAMETER(nout=3) !nout=1 day/nout=2 mth/nout=3 NMC
c
c variables Output
      REAL tnondef(klon,klevSTD,nout)
      REAL tsumSTD(klon,klevSTD,nout)
c
c initialisation variables journalieres en debut de journee
c
      IF(MOD(itap,ecrit_day).EQ.1.) THEN
       DO k=1, nlevSTD
        DO i=1, klon
         tnondef(i,k,1)=0.
         tsumSTD(i,k,1)=0.
        ENDDO !i
       ENDDO !k
      ENDIF
c
c initialisation variables mensuelles (calculees a chaque pas de temps) 
c en debut de mois : nout=2
c
      IF(MOD(itap,ecrit_mth).EQ.1.) THEN
c
       DO k=1, nlevSTD
        DO i=1, klon
         tnondef(i,k,2)=0.
         tsumSTD(i,k,2)=0.
        ENDDO !i
       ENDDO !k
c
c initialisation variables mensuelles - runs type Amip - (calculees toutes les 6h)
c en debut de mois : nout = 3
c
       DO k=1, nlevSTD
        DO i=1, klon
         tnondef(i,k,3)=0.
         tsumSTD(i,k,3)=0.
        ENDDO !i
       ENDDO !k
c
      ENDIF
c
      RETURN
      END  
