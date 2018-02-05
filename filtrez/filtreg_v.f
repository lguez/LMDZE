module filtreg_v_m

  IMPLICIT NONE

contains

  SUBROUTINE filtreg_v(champ, intensive)

    ! From filtrez/filtreg.F, version 1.1.1.1, 2004/05/19 12:53:09
    ! Author: P. Le Van

    ! Matrix filter on longitudes. Matrices have already been
    ! computed. On v grid, direct filter only.

    USE dimens_m, ONLY: iim, jjm
    use filtreg_hemisph_m, only: filtreg_hemisph
    USE inifgn_m, ONLY: sddu, unsddu
    use inifilr_m, only: jfiltnv, jfiltsv, matricevn, matricevs
    use nr_util, only: assert

    REAL, intent(inout):: champ(:, :, :) ! (iim + 1, jjm, :)
    ! en entr\'ee : champ \`a filtrer, en sortie : champ filtr\'e

    logical, intent(in):: intensive
    ! false means the field is weighted by the area of the mesh

    ! Local:
    REAL sdd(iim)

    !-----------------------------------------------------------

    call assert(size(champ, 1) == iim + 1, "filtreg_v iim + 1")
    call assert(size(champ, 2) == jjm, "filtreg_v jjm")

    sdd = merge(sddu, unsddu, intensive)

    call filtreg_hemisph(champ(:, :jfiltnv, :), sdd, matricevn)
    call filtreg_hemisph(champ(:, jfiltsv:jjm, :), sdd, matricevs)

  END SUBROUTINE filtreg_v

end module filtreg_v_m
