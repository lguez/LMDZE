module massbar_m

  IMPLICIT NONE

contains

  SUBROUTINE massbar(masse, massebx, masseby)

    ! From LMDZ4/libf/dyn3d/massbar.F, version 1.1.1.1, 2004/05/19 12:53:05
    ! Authors: P. Le Van, F. Hourdin.

    ! Calcule les moyennes en x et y de la masse d'air dans chaque
    ! maille. Cf. "inigeom.txt" et "massbar.txt".

    USE dimens_m, ONLY: iim, llm
    USE paramet_m, ONLY: iip1, ip1jm, ip1jmp1
    USE comgeom, ONLY: alpha1p2, alpha1p4, alpha2p3, alpha3p4

    REAL, intent(in):: masse(ip1jmp1,llm)
    real, intent(out):: massebx(ip1jmp1,llm), masseby(ip1jm,llm)

    ! Local:
    integer l, ij

    !--------------------------------------------------------------

    DO l = 1, llm
       DO ij = 1, ip1jmp1 - 1
          massebx(ij,l) = masse(ij, l) * alpha1p2(ij) + &
               masse(ij+1, l) * alpha3p4(ij+1)
       ENDDO

       ! correction pour massebx(iip1,j) 
       ! massebx(iip1,j)= massebx(1,j) 
       DO ij = iip1, ip1jmp1, iip1
          massebx(ij,l) = massebx(ij - iim,l)
       ENDDO

       DO ij = 1,ip1jm
          masseby(ij,l) = masse(ij, l) * alpha2p3(ij) + &
               masse(ij+iip1, l) * alpha1p4(ij+iip1)
       ENDDO
    end DO

  END SUBROUTINE massbar

end module massbar_m
