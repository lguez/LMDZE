module covcont_m

  IMPLICIT NONE

contains

  SUBROUTINE covcont(klevel, ucov, vcov, ucont, vcont)

    ! From LMDZ4/libf/dyn3d/covcont.F, version 1.1.1.1 2004/05/19 12:53:07

    ! Auteur:  P. Le Van

    ! Objet: calcul des composantes contravariantes \`a partir des
    ! composantes covariantes

    USE paramet_m, only: ip1jmp1, ip1jm, iip2
    USE comgeom, only: unscu2, unscv2

    INTEGER klevel
    REAL, INTENT (IN) :: ucov(ip1jmp1, klevel), vcov(ip1jm, klevel)
    REAL ucont(ip1jmp1, klevel), vcont(ip1jm, klevel)

    ! Local:
    INTEGER l, ij

    !-------------------------------------------------------------------

    DO l = 1, klevel
       DO ij = iip2, ip1jm
          ucont(ij, l) = ucov(ij, l) * unscu2(ij)
       END DO

       DO ij = 1, ip1jm
          vcont(ij, l) = vcov(ij, l) * unscv2(ij)
       END DO
    END DO

  END SUBROUTINE covcont

end module covcont_m
