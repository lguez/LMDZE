SUBROUTINE filtreg(champ, nlat, nbniv, ifiltre, iaire, griscal, iter)

  ! From filtrez/filtreg.F, version 1.1.1.1 2004/05/19 12:53:09

  ! Auteur: P. Le Van 07/10/97 
  ! Objet: filtre matriciel longitudinal, avec les matrices précalculées
  ! pour l'operateur filtre. 

  USE dimens_m
  USE paramet_m
  USE parafilt
  USE coefils

  IMPLICIT NONE

  INTEGER, intent(in):: nlat ! nombre de latitudes a filtrer
  integer, intent(in):: nbniv ! nombre de niveaux verticaux a filtrer 

  REAL, intent(inout):: champ(iip1, nlat, nbniv)
  ! en entree : champ a filtrer, en sortie : champ filtre

  integer, intent(in):: ifiltre
  !  +1 Transformee directe 
  ! -1 Transformee inverse 
  ! +2 Filtre directe 
  ! -2 Filtre inverse 
  ! Variable Intensive 
  ! ifiltre = 1 filtre directe 
  ! ifiltre =-1 filtre inverse 
  ! Variable Extensive 
  ! ifiltre = 2 filtre directe 
  ! ifiltre =-2 filtre inverse

  integer, intent(in):: iaire
  !  1 si champ intensif 
  ! 2 si champ extensif (pondere par les aires)

  integer, intent(in):: iter
  !  1 filtre simple 

  LOGICAL, intent(in):: griscal

  ! Variables local to the procedure:

  INTEGER jdfil1, jdfil2, jffil1, jffil2, jdfil, jffil
  INTEGER i, j, l, k
  REAL matriceun, matriceus, matricevn, matricevs, matrinvn, matrinvs
  COMMON /matrfil/matriceun(iim, iim, nfilun), matriceus(iim, iim, nfilus), &
       matricevn(iim, iim, nfilvn), matricevs(iim, iim, nfilvs), &
       matrinvn(iim, iim, nfilun), matrinvs(iim, iim, nfilus)
  REAL eignq(iim), sdd1(iim), sdd2(iim)
  INTEGER hemisph

  !-----------------------------------------------------------

  IF (ifiltre==1 .OR. ifiltre==-1) STOP &
       'Pas de transformee simple dans cette version'

  IF (iter==2) THEN
     PRINT *, ' Pas d iteration du filtre dans cette version !', &
          ' Utiliser old_filtreg et repasser !'
     STOP
  END IF

  IF (ifiltre==-2 .AND. .NOT. griscal) THEN
     PRINT *, ' Cette routine ne calcule le filtre inverse que ', &
          ' sur la grille des scalaires !'
     STOP
  END IF

  IF (ifiltre/=2 .AND. ifiltre/=-2) THEN
     PRINT *, ' Probleme dans filtreg car ifiltre NE 2 et NE -2', &
          ' corriger et repasser !'
     STOP
  END IF

  IF (griscal) THEN
     IF (nlat/=jjp1) THEN
        PRINT 1111
        STOP
     ELSE

        IF (iaire==1) THEN
           CALL scopy(iim, sddv, 1, sdd1, 1)
           CALL scopy(iim, unsddv, 1, sdd2, 1)
        ELSE
           CALL scopy(iim, unsddv, 1, sdd1, 1)
           CALL scopy(iim, sddv, 1, sdd2, 1)
        END IF

        jdfil1 = 2
        jffil1 = jfiltnu
        jdfil2 = jfiltsu
        jffil2 = jjm
     END IF
  ELSE
     IF (nlat/=jjm) THEN
        PRINT 2222
        STOP
     ELSE

        IF (iaire==1) THEN
           CALL scopy(iim, sddu, 1, sdd1, 1)
           CALL scopy(iim, unsddu, 1, sdd2, 1)
        ELSE
           CALL scopy(iim, unsddu, 1, sdd1, 1)
           CALL scopy(iim, sddu, 1, sdd2, 1)
        END IF

        jdfil1 = 1
        jffil1 = jfiltnv
        jdfil2 = jfiltsv
        jffil2 = jjm
     END IF
  END IF


  DO hemisph = 1, 2

     IF (hemisph==1) THEN
        jdfil = jdfil1
        jffil = jffil1
     ELSE
        jdfil = jdfil2
        jffil = jffil2
     END IF


     DO l = 1, nbniv
        DO j = jdfil, jffil


           DO i = 1, iim
              champ(i, j, l) = champ(i, j, l)*sdd1(i)
           END DO


           IF (hemisph==1) THEN

              IF (ifiltre==-2) THEN
                 DO k = 1, iim
                    eignq(k) = 0.0
                 END DO
                 DO k = 1, iim
                    DO i = 1, iim
                       eignq(k) = eignq(k) + matrinvn(k, i, j)*champ(i, j, l)
                    END DO
                 END DO
              ELSE IF (griscal) THEN
                 DO k = 1, iim
                    eignq(k) = 0.0
                 END DO
                 DO i = 1, iim
                    DO k = 1, iim
                       eignq(k) = eignq(k) + matriceun(k, i, j)*champ(i, j, l)
                    END DO
                 END DO
              ELSE
                 DO k = 1, iim
                    eignq(k) = 0.0
                 END DO
                 DO i = 1, iim
                    DO k = 1, iim
                       eignq(k) = eignq(k) + matricevn(k, i, j)*champ(i, j, l)
                    END DO
                 END DO
              END IF

           ELSE

              IF (ifiltre==-2) THEN
                 DO k = 1, iim
                    eignq(k) = 0.0
                 END DO
                 DO i = 1, iim
                    DO k = 1, iim
                       eignq(k) = eignq(k) + matrinvs(k, i, j-jfiltsu+1)*champ(i, j, &
                            l)
                    END DO
                 END DO
              ELSE IF (griscal) THEN
                 DO k = 1, iim
                    eignq(k) = 0.0
                 END DO
                 DO i = 1, iim
                    DO k = 1, iim
                       eignq(k) = eignq(k) + matriceus(k, i, j-jfiltsu+1)*champ(i, j &
                            , l)
                    END DO
                 END DO
              ELSE
                 DO k = 1, iim
                    eignq(k) = 0.0
                 END DO
                 DO i = 1, iim
                    DO k = 1, iim
                       eignq(k) = eignq(k) + matricevs(k, i, j-jfiltsv+1)*champ(i, j &
                            , l)
                    END DO
                 END DO
              END IF

           END IF

           IF (ifiltre==2) THEN
              DO i = 1, iim
                 champ(i, j, l) = (champ(i, j, l)+eignq(i))*sdd2(i)
              end DO
           ELSE
              DO i = 1, iim
                 champ(i, j, l) = (champ(i, j, l)-eignq(i))*sdd2(i)
              end DO
           END IF

           champ(iip1, j, l) = champ(1, j, l)

        END DO

     END DO

  end DO

1111 FORMAT (//20X, 'ERREUR dans le dimensionnement du tableau &
       & CHAMP a filtrer, sur la grille des scalaires'/)
2222 FORMAT (//20X, 'ERREUR dans le dimensionnement du tableau &
       &CHAMP a filtrer, sur la grille de V ou de Z'/)

END SUBROUTINE filtreg
