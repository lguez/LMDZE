module grad_m

  IMPLICIT NONE

contains

  SUBROUTINE grad(klevel, pg, pgx, pgy)

    ! From LMDZ4/libf/dyn3d/grad.F, version 1.1.1.1 2004/05/19 12:53:05
    ! P. Le Van

    ! Calcul des composantes covariantes en x et y du gradient de g.

    USE dimens_m, ONLY : iim
    USE paramet_m, ONLY : iip1, ip1jm, ip1jmp1

    INTEGER, intent(in):: klevel
    REAL, intent(in):: pg(ip1jmp1, klevel)
    REAL, intent(out):: pgx(ip1jmp1, klevel) , pgy(ip1jm, klevel)

    ! Local:
    INTEGER l, ij

    !----------------------------------------------------------------

    DO l = 1, klevel
       DO ij = 1, ip1jmp1 - 1
          pgx(ij, l) = pg(ij +1, l) - pg(ij, l)
       end DO

       ! correction pour pgx(ip1, j, l)
       ! pgx(iip1, j, l)= pgx(1, j, l)
       DO ij = iip1, ip1jmp1, iip1
          pgx(ij, l) = pgx(ij -iim, l)
       end DO

       DO ij = 1, ip1jm
          pgy(ij, l) = pg(ij, l) - pg(ij +iip1, l)
       end DO
    end DO

  END SUBROUTINE grad

end module grad_m
