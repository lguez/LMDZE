module laplacien_m

  IMPLICIT NONE

contains

  SUBROUTINE laplacien(klevel, teta)

    ! From LMDZ4/libf/dyn3d/laplacien.F, version 1.1.1.1 2004/05/19 12:53:06
    ! P. Le Van
    ! Calcul de div(grad) de teta.

    use divergf_m, only: divergf
    use filtreg_m, only: filtreg
    use grad_m, only: grad
    USE paramet_m, ONLY: ip1jm, ip1jmp1, jjp1

    INTEGER, intent(in):: klevel
    REAL, intent(inout):: teta(ip1jmp1, klevel)

    ! Variables locales:
    REAL ghy(ip1jm, klevel), ghx(ip1jmp1, klevel)

    !-----------------------------------------------------------------

    CALL filtreg(teta, jjp1, klevel, 2, 1, .TRUE.)
    CALL grad(klevel, teta, ghx, ghy)
    CALL divergf(klevel, ghx, ghy, teta)

  END SUBROUTINE laplacien

end module laplacien_m
