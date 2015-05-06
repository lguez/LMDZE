module gradiv2_m

  IMPLICIT NONE

contains

  SUBROUTINE gradiv2(xcov, ycov, ld, gdx, gdy, cdivu)

    ! From LMDZ4/libf/dyn3d/gradiv2.F, version 1.1.1.1 2004/05/19 12:53:07
    ! P. Le Van
    ! Calcul du gradient de la divergence du vecteur v.

    USE comgeom, ONLY: cuvscvgam1, cvuscugam1, unsair_gam1, unsapolnga1, &
         unsapolsga1
    USE dimens_m, ONLY: iim, jjm
    use divergf_m, only: divergf
    USE filtreg_scal_m, ONLY: filtreg_scal
    use grad_m, only: grad
    use laplacien_m, only: laplacien
    use nr_util, only: assert_eq, assert

    ! Composantes covariantes de v :
    REAL, intent(in):: xcov(:, :, :) ! (iim + 1, jjm + 1, klevel)
    REAL, intent(in):: ycov(:, :, :) ! (iim + 1, jjm, klevel)

    integer, intent(in):: ld
    REAL, intent(out):: gdx(:, :, :) ! (iim + 1, jjm + 1, klevel)
    REAL, intent(out):: gdy(:, :, :) ! (iim + 1, jjm, klevel)
    real, intent(in):: cdivu

    ! Variables locales :
    REAL nugrads, div(iim + 1, jjm + 1, size(xcov, 3))
    INTEGER iter, klevel

    !--------------------------------------------------------------

    call assert((/size(xcov, 1), size(ycov, 1), size(gdx, 1), size(gdy, 1)/) &
         == iim + 1, "gradiv2 iim")
    call assert((/size(xcov, 2) - 1, size(ycov, 2), size(gdx, 2) - 1, &
         size(gdy, 2)/) == jjm, "gradiv2 iim")
    klevel = assert_eq(size(xcov, 3), size(ycov, 3), size(gdx, 3), &
         size(gdy, 3), "gradiv2 klevel")

    CALL divergf(klevel, xcov, ycov, div)

    IF (ld > 1) THEN
       CALL laplacien(klevel, div)

       ! Itération de l'opérateur laplacien_gam
       DO iter = 1, ld - 2
          CALL laplacien_gam(klevel, cuvscvgam1, cvuscugam1, unsair_gam1, &
               unsapolnga1, unsapolsga1, div, div)
       END DO
    ENDIF

    CALL filtreg_scal(div, direct = .true., intensive = .true.)
    CALL grad(klevel, div, gdx, gdy)
    nugrads = (-1.)**ld * cdivu

    gdx = gdx * nugrads
    gdy = gdy * nugrads

  END SUBROUTINE gradiv2

end module gradiv2_m
