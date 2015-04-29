
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/laplacien_gam.F,v 1.1.1.1
! 2004/05/19 12:53:06 lmdzadmin Exp $

SUBROUTINE laplacien_gam(klevel, cuvsga, cvusga, unsaigam, unsapolnga, &
    unsapolsga, teta, divgra)

  ! P. Le Van

  ! ************************************************************

  ! ....   calcul de  (div( grad ))   de   teta  .....
  ! ************************************************************
  ! klevel et teta  sont des arguments  d'entree pour le s-prog
  ! divgra     est  un argument  de sortie pour le s-prog

  USE grad_m, ONLY: grad
  USE dimens_m
  USE paramet_m
  USE comgeom
  IMPLICIT NONE



  ! ............     variables  en arguments    ..........

  INTEGER, INTENT (IN) :: klevel
  REAL teta(ip1jmp1, klevel), divgra(ip1jmp1, klevel)
  REAL cuvsga(ip1jm), cvusga(ip1jmp1), unsaigam(ip1jmp1), unsapolnga, &
    unsapolsga

  ! ...........    variables  locales    .................

  REAL ghy(ip1jm, klevel), ghx(ip1jmp1, klevel)
  ! ......................................................

  CALL scopy(ip1jmp1*klevel, teta, 1, divgra, 1)

  CALL grad(klevel, divgra, ghx, ghy)

  CALL diverg_gam(klevel, cuvsga, cvusga, unsaigam, unsapolnga, unsapolsga, &
    ghx, ghy, divgra)




  RETURN
END SUBROUTINE laplacien_gam
