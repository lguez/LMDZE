module nxgraro2_m

  IMPLICIT NONE

contains

  SUBROUTINE nxgraro2(xcov, ycov, lr, grx, gry, crot)

    ! From LMDZ4/libf/dyn3d/nxgraro2.F, version 1.1.1.1 2004/05/19 12:53:06

    ! P. Le Van
    ! Calcul de nxgrad(rot) du vecteur v

    USE dimens_m, ONLY: iim, jjm
    USE filtreg_v_m, ONLY: filtreg_v
    use nr_util, only: assert, assert_eq
    use rotatf_m, only: rotatf

    ! Composantes covariantes de v :
    REAL, intent(in):: xcov(:, :, :) ! (iim + 1, jjm + 1, :)
    REAL, intent(in):: ycov(:, :, :) ! (iim + 1, jjm, :)

    integer, intent(in):: lr
    REAL, intent(out):: grx(:, :, :) ! (iim + 1, jjm + 1, :)
    REAL, intent(out):: gry(:, :, :) ! (iim + 1, jjm, :)
    real, intent(in):: crot

    ! Variables locales

    INTEGER klevel, iter
    REAL rot(iim + 1, jjm, size(xcov, 3)) , nugradrs

    !----------------------------------------------------------

    call assert((/size(xcov, 1), size(ycov, 1), size(grx, 1), size(gry, 1)/) &
         == iim + 1, "nxgraro2 iim")
    call assert((/size(xcov, 2) - 1, size(ycov, 2), size(grx, 2) - 1, &
         size(gry, 2)/) == jjm, "nxgraro2 jjm")
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
