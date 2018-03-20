module divergf_m

  IMPLICIT NONE

contains

  SUBROUTINE divergf(klevel, x, y, div)

    ! From libf/dyn3d/divergf.F, version 1.1.1.1 2004/05/19 12:53:05

    ! P. Le Van

    ! Calcule la divergence à tous les niveaux d'un vecteur de
    ! composantes x et y. x et y sont des composantes covariantes.

    USE comgeom, ONLY: apoln, apols, cuvsurcv_2d, cvusurcu_2d, unsaire_2d
    USE dimensions, ONLY: iim, jjm
    USE filtreg_scal_m, ONLY: filtreg_scal

    INTEGER, intent(in):: klevel
    REAL, intent(in):: x(iim + 1, jjm + 1, klevel), y(iim + 1, jjm, klevel)
    real, intent(out):: div(iim + 1, jjm + 1, klevel) ! in (unit of x, y) m-2

    ! Variables locales :

    INTEGER l, i, j

    !------------------------------------------------------------

    DO l = 1, klevel
       forall (i = 2:iim + 1, j = 2:jjm) div(i, j, l) = cvusurcu_2d(i, j) &
            * x(i, j, l) - cvusurcu_2d(i - 1, j) * x(i - 1, j , l) &
            + cuvsurcv_2d(i, j - 1) * y(i, j - 1, l) - cuvsurcv_2d(i, j) &
            * y(i, j, l)

       div(1, 2:jjm, l) = div(iim + 1, 2:jjm, l)

       ! Calcul aux pôles 
       div(:, 1, l) = - SUM(cuvsurcv_2d(:iim, 1) * y(:iim, 1, l)) / apoln
       div(:, jjm + 1, l) = SUM(cuvsurcv_2d(:iim, jjm) * y(:iim, jjm, l)) &
            / apols
    end DO

    CALL filtreg_scal(div, direct = .true., intensive = .false.)

    DO l = 1, klevel
       div(:, 2:jjm, l) = div(:, 2:jjm, l) * unsaire_2d(:, 2:jjm)
    ENDDO

  END SUBROUTINE divergf

end module divergf_m
