module massdair_m

  IMPLICIT NONE

contains

  SUBROUTINE massdair(p, masse)

    ! From LMDZ4/libf/dyn3d/massdair.F, version 1.1.1.1 2004/05/19 12:53:07

    ! Calcule la masse d'air dans chaque maille.
    ! Auteurs : P. Le Van , F. Hourdin.

    ! Méthode pour calculer massebx et masseby. A chaque point
    ! scalaire P(i, j) sont affectés quatre coefficients d'aires.

    ! alpha1(i, j) calculé au point (i + 1/4, j - 1/4)
    ! alpha2(i, j) calculé au point (i + 1/4, j + 1/4)
    ! alpha3(i, j) calculé au point (i - 1/4, j + 1/4)
    ! alpha4(i, j) calculé au point (i - 1/4, j - 1/4)

    ! Avec alpha1(i, j) = aire(i + 1/4, j - 1/4)/ aire(i, j) 

    ! Pour plus de détails, voir sous-programme "iniconst" et
    ! "massdair.txt".

    USE comgeom, ONLY: airesurg
    USE dimens_m, ONLY: iim, llm
    USE paramet_m, ONLY: iip1, ip1jmp1, llmp1

    REAL, intent(in):: p(ip1jmp1, llmp1) ! aux interfaces des llm couches
    real, intent(out):: masse(ip1jmp1, llm)

    ! Variables locales
    INTEGER l, ij

    !----------------------------------------------------------

    DO l = 1 , llm
       masse(:, l) = airesurg(:) * (p(:, l) - p(:, l + 1))

       DO ij = 1, ip1jmp1, iip1
          masse(ij + iim, l) = masse(ij, l)
       ENDDO
    end DO

  END SUBROUTINE massdair

end module massdair_m
