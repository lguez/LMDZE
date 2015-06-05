module filtreg_scal_m

  IMPLICIT NONE

contains

  SUBROUTINE filtreg_scal(champ, direct, intensive)

    ! From filtrez/filtreg.F, version 1.1.1.1, 2004/05/19 12:53:09
    ! Author: P. Le Van
    ! Objet : filtre matriciel longitudinal, avec les matrices pr\'ecalcul\'ees
    ! pour l'op\'erateur filtre. 

    USE dimens_m, ONLY: iim, jjm
    use filtreg_hemisph_m, only: filtreg_hemisph
    USE inifgn_m, ONLY: sddv, unsddv
    use inifilr_m, only: jfiltnu, jfiltsu, matriceun, matriceus, matrinvn, &
         matrinvs
    use nr_util, only: assert

    REAL, intent(inout):: champ(:, :, :) ! (iim + 1, jjm + 1, :)
    ! en entr\'ee : champ \`a filtrer, en sortie : champ filtr\'e

    logical, intent(in):: direct ! filtre direct ou inverse 

    logical, intent(in):: intensive
    ! champ intensif ou extensif (pond\'er\'e par les aires)

    ! Local:
    REAL sdd1(iim), sdd2(iim)

    !-----------------------------------------------------------

    call assert(size(champ, 1) == iim + 1, "filtreg_scal iim + 1")
    call assert(size(champ, 2) == jjm + 1, "filtreg_scal jjm + 1")

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

  END SUBROUTINE filtreg_scal

end module filtreg_scal_m
