module gradiv2_m

  IMPLICIT NONE

contains

  SUBROUTINE gradiv2(klevel, xcov, ycov, ld, gdx, gdy, cdivu)

    ! From LMDZ4/libf/dyn3d/gradiv2.F, version 1.1.1.1 2004/05/19 12:53:07
    ! P. Le Van
    ! Calcul de grad div du vecteur v.

    USE dimens_m, ONLY : llm
    USE paramet_m, ONLY : ip1jm, ip1jmp1, jjp1
    USE comgeom, ONLY : cuvscvgam1, cvuscugam1, unsair_gam1, unsapolnga1, &
         unsapolsga1
    USE filtreg_m, ONLY : filtreg

    INTEGER, intent(in):: klevel

    ! composantes covariantes de v:
    REAL, intent(in):: xcov(ip1jmp1, klevel), ycov(ip1jm, klevel)

    integer, intent(in):: ld
    REAL, intent(out):: gdx(ip1jmp1, klevel), gdy(ip1jm, klevel)
    real, intent(in):: cdivu

    ! Variables locales :
    REAL div(ip1jmp1, llm)
    REAL nugrads
    INTEGER iter

    !--------------------------------------------------------------

    gdx = xcov
    gdy = ycov

    CALL divergf(klevel, gdx, gdy, div)

    IF (ld > 1) THEN
       CALL laplacien(klevel, div, div)

       ! Itération de l'opérateur laplacien_gam
       DO iter = 1, ld -2
          CALL laplacien_gam(klevel, cuvscvgam1, cvuscugam1, unsair_gam1, &
               unsapolnga1, unsapolsga1, div, div)
       END DO
    ENDIF

    CALL filtreg(div, jjp1, klevel, 2, 1, .TRUE., 1)
    CALL grad(klevel, div, gdx, gdy)
    nugrads = (-1.)**ld * cdivu

    gdx = gdx * nugrads
    gdy = gdy * nugrads

  END SUBROUTINE gradiv2

end module gradiv2_m
