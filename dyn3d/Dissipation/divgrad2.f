module divgrad2_m

  IMPLICIT NONE

contains

  SUBROUTINE divgrad2(klevel, h, deltapres, lh, divgra, cdivh)

    ! From LMDZ4/libf/dyn3d/divgrad2.F, version 1.1.1.1 2004/05/19 12:53:06
    ! P. Le Van

    ! Calcul de div(grad) de (pext * h)

    USE comgeom, ONLY: cuvscvgam2, cvuscugam2, unsair_gam2, unsapolnga2, &
         unsapolsga2
    use laplacien_gam_m, only: laplacien_gam
    USE laplacien_m, ONLY: laplacien
    USE paramet_m, ONLY: ip1jmp1

    INTEGER, intent(in):: klevel
    REAL, intent(in):: h(ip1jmp1, klevel), deltapres(ip1jmp1, klevel)
    integer, intent(in):: lh
    REAL, intent(out):: divgra(ip1jmp1, klevel)
    real, intent(in):: cdivh

    ! Variables locales
    REAL sqrtps(ip1jmp1, klevel)
    INTEGER iter

    !-----------------------------------------------------------------

    divgra = h
    CALL laplacien(klevel, divgra)
    sqrtps = SQRT(deltapres)
    divgra = divgra * sqrtps

    ! Itération de l'opérateur laplacien_gam
    DO iter = 1, lh - 2
       CALL laplacien_gam(klevel, cuvscvgam2, cvuscugam2, unsair_gam2, &
            unsapolnga2, unsapolsga2, divgra, divgra)
    ENDDO

    divgra = divgra * sqrtps
    CALL laplacien(klevel, divgra)
    divgra = (-1.)**lh * cdivh * divgra / deltapres

  END SUBROUTINE divgrad2

end module divgrad2_m
