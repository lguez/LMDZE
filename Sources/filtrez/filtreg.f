module filtreg_m

  IMPLICIT NONE

contains

  SUBROUTINE filtreg(champ, direct, intensive)

    ! From filtrez/filtreg.F, version 1.1.1.1, 2004/05/19 12:53:09
    ! Author: P. Le Van
    ! Objet : filtre matriciel longitudinal, avec les matrices pr\'ecalcul\'ees
    ! pour l'op\'erateur filtre. 

    USE coefils, ONLY: sddu, sddv, unsddu, unsddv
    USE dimens_m, ONLY: iim, jjm
    use inifilr_m, only: jfiltnu, jfiltnv, jfiltsu, jfiltsv, matriceun, &
         matriceus, matricevn, matricevs, matrinvn, matrinvs
    use nr_util, only: assert

    REAL, intent(inout):: champ(:, :, :) ! (iim + 1, nlat, nbniv)
    ! en entr\'ee : champ \`a filtrer, en sortie : champ filtr\'e

    logical, intent(in):: direct ! filtre direct ou inverse 

    logical, intent(in):: intensive
    ! champ intensif ou extensif (pond\'er\'e par les aires)

    ! Local:
    LOGICAL griscal
    INTEGER nlat ! nombre de latitudes \`a filtrer
    integer nbniv ! nombre de niveaux verticaux \`a filtrer 
    INTEGER jdfil1, jdfil2, jffil1, jffil2, jdfil, jffil
    INTEGER j, l
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

       DO l = 1, nbniv
          DO j = jdfil, jffil
             champ(:iim, j, l) = champ(:iim, j, l) * sdd1

             IF (hemisph==1) THEN
                IF (.not. direct) THEN
                   eignq = matmul(matrinvn(:, :, j), champ(:iim, j, l))
                ELSE IF (griscal) THEN
                   eignq = matmul(matriceun(:, :, j), champ(:iim, j, l))
                ELSE
                   eignq = matmul(matricevn(:, :, j), champ(:iim, j, l))
                END IF
             ELSE
                IF (.not. direct) THEN
                   eignq = matmul(matrinvs(:, :, j - jfiltsu + 1), &
                        champ(:iim, j, l))
                ELSE IF (griscal) THEN
                   eignq = matmul(matriceus(:, :, j - jfiltsu + 1), &
                        champ(:iim, j, l))
                ELSE
                   eignq = matmul(matricevs(:, :, j - jfiltsv + 1), &
                        champ(:iim, j, l))
                END IF
             END IF

             IF (direct) THEN
                champ(:iim, j, l) = (champ(:iim, j, l) + eignq) * sdd2
             ELSE
                champ(:iim, j, l) = (champ(:iim, j, l) - eignq) * sdd2
             END IF

             champ(iim + 1, j, l) = champ(1, j, l)
          END DO
       END DO
    end DO

  END SUBROUTINE filtreg

end module filtreg_m
