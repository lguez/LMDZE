module massbarxy_m

  IMPLICIT NONE

contains

  SUBROUTINE massbarxy(masse, massebxy)

    ! From LMDZ4/libf/dyn3d/massbarxy.F, version 1.1.1.1 2004/05/19 12:53:07

    ! Calcule les moyennes en x et y de la masse d'air dans chaque maille.
    ! Authors: P. Le Van, F. Hourdin.

    USE comgeom, ONLY: alpha1, alpha2, alpha3, alpha4
    USE dimens_m, ONLY: iim, llm
    USE paramet_m, ONLY: iip1, iip2, ip1jm, ip1jmp1

    REAL, intent(in):: masse(ip1jmp1, llm)
    real, intent(out):: massebxy(ip1jm, llm)

    ! Local:
    integer l, ij

    !-------------------------------------------------------------------

    DO l = 1, llm
       DO ij = 1, ip1jm - 1
          massebxy(ij, l) = masse(ij, l) * alpha2(ij) &
               + masse(ij + 1, l) * alpha3(ij + 1) &
               + masse(ij + iip1, l) * alpha1(ij + iip1) &
               + masse(ij + iip2, l) * alpha4(ij + iip2)
       end DO

       !  correction pour massebxy(iip1, j) 
       DO ij = iip1, ip1jm, iip1
          massebxy(ij, l) = massebxy(ij - iim, l)
       end DO
    end DO

  END SUBROUTINE massbarxy

end module massbarxy_m
