module filtreg_m

  IMPLICIT NONE

contains

  SUBROUTINE filtreg(champ, direct, intensive)

    ! From filtrez/filtreg.F, version 1.1.1.1, 2004/05/19 12:53:09
    ! Author: P. Le Van
    ! Objet : filtre matriciel longitudinal, avec les matrices précalculées
    ! pour l'opérateur filtre. 

    USE coefils, ONLY: sddu, sddv, unsddu, unsddv
    USE dimens_m, ONLY: iim, jjm
    use inifilr_m, only: jfiltnu, jfiltnv, jfiltsu, jfiltsv, matriceun, &
         matriceus, matricevn, matricevs, matrinvn, matrinvs
    use nr_util, only: assert

    REAL, intent(inout):: champ(:, :, :) ! (iim + 1, nlat, nbniv)
    ! en entrée : champ à filtrer, en sortie : champ filtré

    logical, intent(in):: direct ! filtre direct ou inverse 

    logical, intent(in):: intensive
    ! champ intensif ou extensif (pondéré par les aires)

    ! Local:
    LOGICAL griscal
    INTEGER nlat ! nombre de latitudes à filtrer
    integer nbniv ! nombre de niveaux verticaux à filtrer 
    INTEGER jdfil1, jdfil2, jffil1, jffil2, jdfil, jffil
    INTEGER i, j, l, k
    REAL eignq(iim), sdd1(iim), sdd2(iim)
    INTEGER hemisph

    !-----------------------------------------------------------

    call assert(size(champ, 1) == iim + 1, "filtreg iim + 1")
    nlat = size(champ, 2)
    nbniv = size(champ, 3)
    call assert(nlat == jjm .or. nlat == jjm + 1, "filtreg nlat")
    griscal = nlat == jjm + 1

    IF (.not. direct .AND. nlat == jjm) THEN
       PRINT *, 'filtreg: inverse filter on scalar grid only'
       STOP 1
    END IF

    IF (griscal) THEN
       IF (intensive) THEN
          sdd1 = sddv
          sdd2 = unsddv
       ELSE
          sdd1 = unsddv
          sdd2 = sddv
       END IF

       jdfil1 = 2
       jffil1 = jfiltnu
       jdfil2 = jfiltsu
       jffil2 = jjm
    ELSE
       IF (intensive) THEN
          sdd1 = sddu
          sdd2 = unsddu
       ELSE
          sdd1 = unsddu
          sdd2 = sddu
       END IF

       jdfil1 = 1
       jffil1 = jfiltnv
       jdfil2 = jfiltsv
       jffil2 = jjm
    END IF

    DO hemisph = 1, 2
       IF (hemisph==1) THEN
          jdfil = jdfil1
          jffil = jffil1
       ELSE
          jdfil = jdfil2
          jffil = jffil2
       END IF

       loop_vertical: DO l = 1, nbniv
          loop_latitude: DO j = jdfil, jffil
             DO i = 1, iim
                champ(i, j, l) = champ(i, j, l)*sdd1(i)
             END DO

             IF (hemisph==1) THEN
                IF (.not. direct) THEN
                   DO k = 1, iim
                      eignq(k) = 0.
                   END DO
                   DO k = 1, iim
                      DO i = 1, iim
                         eignq(k) = eignq(k) + matrinvn(k, i, j)*champ(i, j, l)
                      END DO
                   END DO
                ELSE IF (griscal) THEN
                   DO k = 1, iim
                      eignq(k) = 0.
                   END DO
                   DO i = 1, iim
                      DO k = 1, iim
                         eignq(k) = eignq(k) + matriceun(k, i, j) &
                              * champ(i, j, l)
                      END DO
                   END DO
                ELSE
                   DO k = 1, iim
                      eignq(k) = 0.
                   END DO
                   DO i = 1, iim
                      DO k = 1, iim
                         eignq(k) = eignq(k) + matricevn(k, i, j) &
                              * champ(i, j, l)
                      END DO
                   END DO
                END IF
             ELSE
                IF (.not. direct) THEN
                   DO k = 1, iim
                      eignq(k) = 0.
                   END DO
                   DO i = 1, iim
                      DO k = 1, iim
                         eignq(k) = eignq(k) + matrinvs(k, i, j-jfiltsu+1) &
                              *champ(i, j, l)
                      END DO
                   END DO
                ELSE IF (griscal) THEN
                   DO k = 1, iim
                      eignq(k) = 0.
                   END DO
                   DO i = 1, iim
                      DO k = 1, iim
                         eignq(k) = eignq(k) + matriceus(k, i, j-jfiltsu+1) &
                              *champ(i, j , l)
                      END DO
                   END DO
                ELSE
                   DO k = 1, iim
                      eignq(k) = 0.
                   END DO
                   DO i = 1, iim
                      DO k = 1, iim
                         eignq(k) = eignq(k) + matricevs(k, i, j-jfiltsv+1) &
                              *champ(i, j , l)
                      END DO
                   END DO
                END IF
             END IF

             IF (direct) THEN
                DO i = 1, iim
                   champ(i, j, l) = (champ(i, j, l)+eignq(i))*sdd2(i)
                end DO
             ELSE
                DO i = 1, iim
                   champ(i, j, l) = (champ(i, j, l)-eignq(i))*sdd2(i)
                end DO
             END IF

             champ(iim + 1, j, l) = champ(1, j, l)
          END DO loop_latitude
       END DO loop_vertical
    end DO

  END SUBROUTINE filtreg

end module filtreg_m
