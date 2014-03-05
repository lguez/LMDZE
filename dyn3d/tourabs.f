module tourabs_m

  IMPLICIT NONE

contains

  SUBROUTINE tourabs(ntetaSTD, vcov, ucov, vorabs)

    ! Author: I. Musat, 28 October 2004, adaptation de tourpot
    ! Objet : calcul de la vorticité absolue

    USE dimens_m, ONLY: iim, jjm
    USE paramet_m, ONLY: iip1, ip1jm, ip1jmp1
    USE comconst, ONLY: rad
    USE comgeom, ONLY: cu, cv, fext, rlatv, unsairez
    USE filtreg_m, ONLY: filtreg
    USE nr_util, ONLY: pi

    INTEGER, intent(in):: ntetaSTD
    REAL, intent(in):: vcov(ip1jm, ntetaSTD), ucov(ip1jmp1, ntetaSTD)

    REAL, intent(out):: vorabs(ip1jm, ntetaSTD)
    ! vorabs = Filtre(d(vcov)/dx - d(ucov)/dy) + fext

    ! Variables locales:
    INTEGER l, ij, i, j
    REAL rot(ip1jm, ntetaSTD)

    !--------------------------------------------------------------------

    ! Calcul du rotationnel du vent puis filtrage

    DO l = 1, ntetaSTD
       DO i = 1, iip1
          DO j = 1, jjm
             ij = i + (j - 1) * iip1
             IF (ij <= ip1jm - 1) THEN
                IF (cv(ij) == 0. .OR. cv(ij+1) == 0. .OR. cu(ij) == 0. &
                     .OR. cu(ij+iip1) == 0.) THEN
                   rot(ij, l) = 0.
                ELSE
                   rot(ij, l) = (vcov(ij + 1, l) / cv(ij + 1) - vcov(ij, l) &
                        / cv(ij)) / (2. * pi * RAD * cos(rlatv(j))) &
                        * real(iim) + (ucov(ij + iip1, l) / cu(ij + iip1) &
                        - ucov(ij, l) / cu(ij)) / (pi * RAD) * (real(jjm) - 1.)
                ENDIF
             ENDIF
          end DO
       end DO

       ! correction pour rot(iip1, j, l) 
       ! rot(iip1, j, l) = rot(1, j, l) 
       DO ij = iip1, ip1jm, iip1
          rot(ij, l) = rot(ij - iim, l)
       end DO
    end DO

    CALL filtreg(rot, jjm, ntetaSTD, 2, 1, .FALSE.)

    DO l = 1, ntetaSTD
       DO ij = 1, ip1jm - 1
          vorabs(ij, l) = (rot(ij, l) + fext(ij) * unsairez(ij))
       end DO

       ! correction pour vorabs(iip1, j, l) 
       ! vorabs(iip1, j, l)= vorabs(1, j, l) 
       DO ij = iip1, ip1jm, iip1
          vorabs(ij, l) = vorabs(ij - iim, l)
       end DO
    end DO

  END SUBROUTINE tourabs

end module tourabs_m
