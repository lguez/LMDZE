module laplacien_gam_m

  IMPLICIT NONE

contains

  SUBROUTINE laplacien_gam(klevel, cuvsga, cvusga, unsaigam, unsapolnga, &
       unsapolsga, teta, divgra)

    ! From LMDZ4/libf/dyn3d/laplacien_gam.F,v 1.1.1.1, 2004/05/19 12:53:06

    ! P. Le Van

    ! ************************************************************

    ! ....   calcul de  (div( grad ))   de   teta  .....
    ! ************************************************************
    ! klevel et teta  sont des arguments  d'entree pour le s-prog
    ! divgra     est  un argument  de sortie pour le s-prog

    USE comgeom
    USE dimensions
    use diverg_gam_m, only: diverg_gam
    USE grad_m, ONLY: grad
    USE paramet_m


    ! ............     variables  en arguments    ..........

    INTEGER, INTENT (IN) :: klevel
    REAL teta(ip1jmp1, klevel), divgra(ip1jmp1, klevel)
    REAL cuvsga(ip1jm), cvusga(ip1jmp1), unsaigam(ip1jmp1), unsapolnga, &
         unsapolsga

    ! ...........    variables  locales    .................

    REAL ghy(ip1jm, klevel), ghx(ip1jmp1, klevel)
    ! ......................................................

    divgra = teta
    CALL grad(klevel, divgra, ghx, ghy)
    CALL diverg_gam(klevel, cuvsga, cvusga, unsaigam, unsapolnga, unsapolsga, &
         ghx, ghy, divgra)

  END SUBROUTINE laplacien_gam

end module laplacien_gam_m
