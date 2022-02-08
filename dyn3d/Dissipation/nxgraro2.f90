module nxgraro2_m

  IMPLICIT NONE

contains

  SUBROUTINE nxgraro2(xcov, ycov, lr, grx, gry, crot)

    ! From LMDZ4/libf/dyn3d/nxgraro2.F, version 1.1.1.1, 2004/05/19 12:53:06

    ! P. Le Van
    ! Calcul de nxgrad(rot) du vecteur (xcov, ycov).

    use jumble, only: assert, assert_eq

    USE filtreg_v_m, ONLY: filtreg_v
    use nxgrad_m, only: nxgrad
    use rotatf_m, only: rotatf

    ! Composantes covariantes :
    REAL, intent(in):: xcov(:, :, :) ! (iim + 1, jjm + 1, klevel)
    REAL, intent(in):: ycov(:, :, :) ! (iim + 1, jjm, klevel)

    integer, intent(in):: lr
    REAL, intent(out):: grx(:, :, :) ! (iim + 1, jjm + 1, klevel)
    REAL, intent(out):: gry(:, :, :) ! (iim + 1, jjm, klevel)
    real, intent(in):: crot

    ! Local:
    INTEGER klevel, iter
    REAL rot(size(ycov, 1), size(ycov, 2), size(ycov, 3)) , nugradrs

    !----------------------------------------------------------

    call assert(size(xcov, 1) == [size(ycov, 1), size(grx, 1), size(gry, 1)], &
         "nxgraro2 iim")
    call assert(size(xcov, 2) - 1 == [size(ycov, 2), size(grx, 2) - 1, &
         size(gry, 2)], "nxgraro2 jjm")
    klevel = assert_eq(size(xcov, 3), size(ycov, 3), size(grx, 3), &
         size(gry, 3), "nxgraro2 klevel")

    grx = xcov
    gry = ycov

    CALL rotatf(klevel, grx, gry, rot)
    CALL laplacien_rot(klevel, rot, rot, grx, gry)

    ! Itération de l'opérateur laplacien_rotgam
    DO iter = 1, lr - 2
       CALL laplacien_rotgam(klevel, rot, rot)
    ENDDO

    CALL filtreg_v(rot, intensive = .true.)
    CALL nxgrad(klevel, rot, grx, gry)

    nugradrs = (-1.)**lr * crot
    grx = grx * nugradrs
    gry = gry * nugradrs

  END SUBROUTINE nxgraro2

end module nxgraro2_m
