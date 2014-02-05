module tourpot_m

  IMPLICIT NONE

contains

  SUBROUTINE tourpot(vcov, ucov, massebxy, vorpot)

    ! From LMDZ4/libf/dyn3d/tourpot.F, version 1.1.1.1 2004/05/19 12:53:06

    ! Auteur : P. Le Van
    ! Objet : calcul du tourbillon potentiel

    USE dimens_m, ONLY: iim, jjm, llm
    USE paramet_m, ONLY: iip1, ip1jm, ip1jmp1
    USE comgeom, ONLY: fext
    use filtreg_m, only: filtreg

    REAL, intent(in):: vcov(ip1jm, llm), ucov(ip1jmp1, llm)
    REAL, intent(in):: massebxy(ip1jm, llm)

    real, intent(out):: vorpot(ip1jm, llm)
    ! = (Filtre(d(vcov)/dx - d(ucov)/dy) + fext) / massebxy

    ! Local:
    REAL rot(ip1jm, llm)
    INTEGER l, ij

    !---------------------------------------------------------------

    ! Calcul du rotationnel du vent puis filtrage

    DO l = 1, llm
       DO ij = 1, ip1jm - 1
          rot(ij, l) = vcov(ij + 1, l) - vcov(ij, l) + ucov(ij + iip1, l) &
               - ucov(ij, l)
       end DO

       ! correction pour rot(iip1, j, l)
       ! rot(iip1, j, l) = rot(1, j, l)
       DO ij = iip1, ip1jm, iip1
          rot(ij, l) = rot(ij - iim, l)
       end DO
    end DO

    CALL filtreg(rot, jjm, llm, 2, 1, .FALSE.)

    DO l = 1, llm
       DO ij = 1, ip1jm - 1
          vorpot(ij, l) = (rot(ij, l) + fext(ij)) / massebxy(ij, l)
       end DO

       ! correction pour vorpot(iip1, j, l)
       ! vorpot(iip1, j, l)= vorpot(1, j, l)
       DO ij = iip1, ip1jm, iip1
          vorpot(ij, l) = vorpot(ij - iim, l)
       end DO
    end DO

  END SUBROUTINE tourpot

end module tourpot_m
