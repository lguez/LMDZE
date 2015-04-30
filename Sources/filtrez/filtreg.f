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
    use filtreg_hemisph_m, only: filtreg_hemisph
    use inifilr_m, only: jfiltnu, jfiltnv, jfiltsu, jfiltsv, matriceun, &
         matriceus, matricevn, matricevs, matrinvn, matrinvs
    use nr_util, only: assert

    REAL, intent(inout):: champ(:, :, :) ! (iim + 1, nlat, :)
    ! en entr\'ee : champ \`a filtrer, en sortie : champ filtr\'e

    logical, intent(in):: direct ! filtre direct ou inverse 

    logical, intent(in):: intensive
    ! champ intensif ou extensif (pond\'er\'e par les aires)

    ! Local:
    INTEGER nlat ! nombre de latitudes \`a filtrer
    REAL sdd1(iim), sdd2(iim)

    !-----------------------------------------------------------

    call assert(size(champ, 1) == iim + 1, "filtreg iim + 1")
    nlat = size(champ, 2)
    call assert(nlat == jjm .or. nlat == jjm + 1, "filtreg nlat")

    if (nlat == jjm + 1) then
       IF (intensive) THEN
          sdd1 = sddv
          sdd2 = unsddv
       ELSE
          sdd1 = unsddv
          sdd2 = sddv
       END IF
       if (direct) then
          call filtreg_hemisph(champ(:, 2:jfiltnu, :), sdd1, sdd2, matriceun)
          call filtreg_hemisph(champ(:, jfiltsu:jjm, :), sdd1, sdd2, matriceus)
       else
          call filtreg_hemisph(champ(:, 2:jfiltnu, :), sdd1, sdd2, - matrinvn)
          call filtreg_hemisph(champ(:, jfiltsu:jjm, :), sdd1, sdd2, - matrinvs)
       end if
    else
       IF (intensive) THEN
          sdd1 = sddu
          sdd2 = unsddu
       ELSE
          sdd1 = unsddu
          sdd2 = sddu
       END IF
       if (direct) then
          call filtreg_hemisph(champ(:, :jfiltnv, :), sdd1, sdd2, matricevn)
          call filtreg_hemisph(champ(:, jfiltsv:jjm, :), sdd1, sdd2, matricevs)
       else
          PRINT *, 'filtreg: inverse filter on scalar grid only'
          STOP 1
       END IF
    end if

  END SUBROUTINE filtreg

end module filtreg_m
