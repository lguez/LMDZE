module filtreg_v_m

  IMPLICIT NONE

contains

  SUBROUTINE filtreg_v(champ, intensive)

    ! From filtrez/filtreg.F, version 1.1.1.1, 2004/05/19 12:53:09
    ! Author: P. Le Van

    ! Filtre matriciel longitudinal, avec les matrices
    ! pr\'ecalcul\'ees pour l'op\'erateur filtre. On v grid, direct
    ! only.

    USE coefils, ONLY: sddu, unsddu
    USE dimens_m, ONLY: iim, jjm
    use filtreg_hemisph_m, only: filtreg_hemisph
    use inifilr_m, only: jfiltnv, jfiltsv, matricevn, matricevs
    use nr_util, only: assert

    REAL, intent(inout):: champ(:, :, :) ! (iim + 1, jjm, :)
    ! en entr\'ee : champ \`a filtrer, en sortie : champ filtr\'e

    logical, intent(in):: intensive
    ! champ intensif ou extensif (pond\'er\'e par les aires)

    ! Local:
    REAL sdd1(iim), sdd2(iim)

    !-----------------------------------------------------------

    call assert(size(champ, 1) == iim + 1, "filtreg_v iim + 1")
    call assert(size(champ, 2) == jjm, "filtreg_v jjm")

    IF (intensive) THEN
       sdd1 = sddu
       sdd2 = unsddu
    ELSE
       sdd1 = unsddu
       sdd2 = sddu
    END IF

    call filtreg_hemisph(champ(:, :jfiltnv, :), sdd1, sdd2, matricevn)
    call filtreg_hemisph(champ(:, jfiltsv:jjm, :), sdd1, sdd2, matricevs)

  END SUBROUTINE filtreg_v

end module filtreg_v_m
