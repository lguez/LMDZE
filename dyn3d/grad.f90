module grad_m

  IMPLICIT NONE

contains

  SUBROUTINE grad(klevel, g, gx, gy)

    ! From LMDZ4/libf/dyn3d/grad.F, version 1.1.1.1 2004/05/19 12:53:05
    ! P. Le Van

    ! Calcul des composantes covariantes en x et y du gradient de g.

    USE dimens_m, ONLY: iim
    USE paramet_m, ONLY: iip1, ip1jm, ip1jmp1

    INTEGER, intent(in):: klevel
    REAL, intent(in):: g(ip1jmp1, klevel)
    REAL, intent(out):: gx(ip1jmp1, klevel) , gy(ip1jm, klevel)

    ! Local:
    INTEGER l, ij

    !----------------------------------------------------------------

    DO l = 1, klevel
       DO ij = 1, ip1jmp1 - 1
          gx(ij, l) = g(ij +1, l) - g(ij, l)
       end DO

       ! correction pour gx(ip1, j, l)
       ! gx(iip1, j, l)= gx(1, j, l)
       DO ij = iip1, ip1jmp1, iip1
          gx(ij, l) = gx(ij -iim, l)
       end DO

       DO ij = 1, ip1jm
          gy(ij, l) = g(ij, l) - g(ij +iip1, l)
       end DO
    end DO

  END SUBROUTINE grad

end module grad_m
