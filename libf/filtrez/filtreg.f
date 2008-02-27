!
! $Header: /home/cvsroot/LMDZ4/libf/filtrez/filtreg.F,v 1.1.1.1 2004/05/19 12:53:09 lmdzadmin Exp $
!
      SUBROUTINE filtreg ( champ, nlat, nbniv, ifiltre,iaire,
     .   griscal ,iter)

      use dimens_m
      use paramet_m
      IMPLICIT NONE
c=======================================================================
c
c   Auteur: P. Le Van        07/10/97
c   ------
c
c   Objet: filtre matriciel longitudinal ,avec les matrices precalculees
c                     pour l'operateur  Filtre    .
c   ------
c
c   Arguments:
c   ----------
c
c      nblat                 nombre de latitudes a filtrer
c      nbniv                 nombre de niveaux verticaux a filtrer
c      champ(iip1,nblat,nbniv)  en entree : champ a filtrer
c                            en sortie : champ filtre
c      ifiltre               +1  Transformee directe
c                            -1  Transformee inverse
c                            +2  Filtre directe
c                            -2  Filtre inverse
c
c      iaire                 1   si champ intensif
c                            2   si champ extensif (pondere par les aires)
c
c      iter                  1   filtre simple
c
c=======================================================================
c
c
c                      Variable Intensive
c                ifiltre = 1     filtre directe
c                ifiltre =-1     filtre inverse
c
c                      Variable Extensive
c                ifiltre = 2     filtre directe
c                ifiltre =-2     filtre inverse
c
c
      include "parafilt.h"
      include "coefils.h"
c
      INTEGER nlat,nbniv,ifiltre,iter
      INTEGER i,j,l,k
      INTEGER iim2,immjm
      INTEGER jdfil1,jdfil2,jffil1,jffil2,jdfil,jffil

      REAL  champ( iip1,nlat,nbniv)
      REAL matriceun,matriceus,matricevn,matricevs,matrinvn,matrinvs
      COMMON/matrfil/matriceun(iim,iim,nfilun),matriceus(iim,iim,nfilus)
     ,             , matricevn(iim,iim,nfilvn),matricevs(iim,iim,nfilvs)
     ,             ,  matrinvn(iim,iim,nfilun),matrinvs (iim,iim,nfilus)
      REAL  eignq(iim), sdd1(iim),sdd2(iim)
      LOGICAL    griscal
      INTEGER    hemisph, iaire
c

      IF(ifiltre.EQ.1.or.ifiltre.EQ.-1) 
     *    STOP'Pas de transformee simple dans cette version'

      IF( iter.EQ. 2 )  THEN
       PRINT *,' Pas d iteration du filtre dans cette version !'
     * , ' Utiliser old_filtreg et repasser !'
           STOP
      ENDIF

      IF( ifiltre.EQ. -2 .AND..NOT.griscal )     THEN
       PRINT *,' Cette routine ne calcule le filtre inverse que ',
     * ' sur la grille des scalaires !'
           STOP
      ENDIF

      IF( ifiltre.NE.2 .AND.ifiltre.NE. - 2 )  THEN
       PRINT *,' Probleme dans filtreg car ifiltre NE 2 et NE -2'
     *,' corriger et repasser !'
           STOP
      ENDIF
c

      iim2   = iim * iim
      immjm  = iim * jjm
c
c
      IF( griscal )   THEN
         IF( nlat. NE. jjp1 )  THEN
             PRINT  1111
             STOP
         ELSE
c
             IF( iaire.EQ.1 )  THEN
                CALL SCOPY(  iim,    sddv, 1,  sdd1, 1 ) 
                CALL SCOPY(  iim,  unsddv, 1,  sdd2, 1 )
             ELSE
                CALL SCOPY(  iim,  unsddv, 1,  sdd1, 1 )
                CALL SCOPY(  iim,    sddv, 1,  sdd2, 1 )
             END IF
c
             jdfil1 = 2
             jffil1 = jfiltnu
             jdfil2 = jfiltsu
             jffil2 = jjm
          END IF
      ELSE
          IF( nlat.NE.jjm )  THEN
             PRINT  2222
             STOP
          ELSE
c
             IF( iaire.EQ.1 )  THEN
                CALL SCOPY(  iim,    sddu, 1,  sdd1, 1 ) 
                CALL SCOPY(  iim,  unsddu, 1,  sdd2, 1 )
             ELSE
                CALL SCOPY(  iim,  unsddu, 1,  sdd1, 1 )
                CALL SCOPY(  iim,    sddu, 1,  sdd2, 1 )
             END IF
c
             jdfil1 = 1
             jffil1 = jfiltnv
             jdfil2 = jfiltsv
             jffil2 = jjm
          END IF
      END IF
c
c
      DO 100  hemisph = 1, 2
c
      IF ( hemisph.EQ.1 )  THEN
          jdfil = jdfil1
          jffil = jffil1
      ELSE
          jdfil = jdfil2
          jffil = jffil2
      END IF

 
      DO 50  l = 1, nbniv
      DO 30  j = jdfil,jffil
 
 
      DO  5  i = 1, iim
      champ(i,j,l) = champ(i,j,l) * sdd1(i)
   5  CONTINUE
c

      IF( hemisph. EQ. 1 )      THEN

        IF( ifiltre. EQ. -2 )   THEN
      DO k = 1, iim
         eignq(k) = 0.0
      ENDDO
      DO k = 1, iim
      DO i = 1, iim
         eignq(k) = eignq(k) + matrinvn(k,i,j)*champ(i,j,l)
      ENDDO
      ENDDO
        ELSE IF ( griscal )     THEN
      DO k = 1, iim
         eignq(k) = 0.0
      ENDDO
      DO i = 1, iim
      DO k = 1, iim
         eignq(k) = eignq(k) + matriceun(k,i,j)*champ(i,j,l)
      ENDDO
      ENDDO
        ELSE 
      DO k = 1, iim
         eignq(k) = 0.0
      ENDDO
      DO i = 1, iim
      DO k = 1, iim
         eignq(k) = eignq(k) + matricevn(k,i,j)*champ(i,j,l)
      ENDDO
      ENDDO
        ENDIF

      ELSE

        IF( ifiltre. EQ. -2 )   THEN
      DO k = 1, iim
         eignq(k) = 0.0
      ENDDO
      DO i = 1, iim
      DO k = 1, iim
         eignq(k) = eignq(k) + matrinvs(k,i,j-jfiltsu+1)*champ(i,j,l)
      ENDDO
      ENDDO
        ELSE IF ( griscal )     THEN
      DO k = 1, iim
         eignq(k) = 0.0
      ENDDO
      DO i = 1, iim
      DO k = 1, iim
         eignq(k) = eignq(k) + matriceus(k,i,j-jfiltsu+1)*champ(i,j,l)
      ENDDO
      ENDDO
        ELSE 
      DO k = 1, iim
         eignq(k) = 0.0
      ENDDO
      DO i = 1, iim
      DO k = 1, iim
         eignq(k) = eignq(k) + matricevs(k,i,j-jfiltsv+1)*champ(i,j,l)
      ENDDO
      ENDDO
        ENDIF

      ENDIF
c
      IF( ifiltre.EQ. 2 )  THEN
        DO 15 i = 1, iim
        champ( i,j,l ) = ( champ(i,j,l) + eignq(i) ) * sdd2(i)
  15    CONTINUE
      ELSE
        DO 16 i=1,iim
        champ( i,j,l ) = ( champ(i,j,l) - eignq(i) ) * sdd2(i)
16      CONTINUE
      ENDIF
c
      champ( iip1,j,l ) = champ( 1,j,l )
c
  30  CONTINUE
c
  50  CONTINUE
c    
 100  CONTINUE
c
1111  FORMAT(//20x,'ERREUR dans le dimensionnement du tableau  CHAMP a 
     *filtrer, sur la grille des scalaires'/)
2222  FORMAT(//20x,'ERREUR dans le dimensionnement du tableau CHAMP a fi
     *ltrer, sur la grille de V ou de Z'/)
      RETURN
      END
