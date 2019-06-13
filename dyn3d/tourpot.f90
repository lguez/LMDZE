module tourpot_m

  IMPLICIT NONE

contains

  SUBROUTINE tourpot(vcov, ucov, massebxy, vorpot)

    ! From LMDZ4/libf/dyn3d/tourpot.F, version 1.1.1.1, 2004/05/19 12:53:06

    ! Author: P. Le Van
    ! Objet : calcul du tourbillon potentiel

    USE comgeom, ONLY: fext_2d
    USE dimensions, ONLY: iim, jjm, llm
    use filtreg_v_m, only: filtreg_v

    REAL, intent(in):: vcov(:, :, :) ! (iim + 1, jjm, llm)
    REAL, intent(in):: ucov(:, :, :) ! (iim + 1, jjm + 1, llm)
    REAL, intent(in):: massebxy(:, :, :) ! (iim + 1, jjm, llm) mass of grid cell

    real, intent(out):: vorpot(:, :, :) ! (iim + 1, jjm, llm)
    ! = (Filtre(d(vcov)/dx - d(ucov)/dy) + fext) / massebxy

    ! Local:

    REAL rot(iim + 1, jjm, llm)
    ! relative vorticity multiplied by cell area, in m2 s-1

    INTEGER l, i, j

    !---------------------------------------------------------------

    ! Calcul du rotationnel du vent puis filtrage

    forall (i = 1: iim, j = 1: jjm) rot(i, j, :) &
         = vcov(i + 1, j, :) - vcov(i, j, :) + ucov(i, j + 1, :) - ucov(i, j, :)
    rot(iim + 1, :, :) = rot(1, :, :)

    CALL filtreg_v(rot, intensive = .true.)

    forall (l = 1: llm) vorpot(:iim, :, l) &
         = (rot(:iim, :, l) + fext_2d(:iim, :)) / massebxy(:iim, :, l)
    vorpot(iim + 1, :, :)= vorpot(1, :, :)

  END SUBROUTINE tourpot

end module tourpot_m
