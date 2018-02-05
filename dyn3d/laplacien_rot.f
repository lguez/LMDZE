
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/laplacien_rot.F,v 1.1.1.1
! 2004/05/19 12:53:07 lmdzadmin Exp $

SUBROUTINE laplacien_rot(klevel, rotin, rotout, ghx, ghy)

  ! P. Le Van

  ! ************************************************************
  ! ...  calcul de  ( rotat x nxgrad )  du rotationnel rotin  .
  ! ************************************************************

  ! klevel et rotin  sont des arguments  d'entree pour le s-prog
  ! rotout           est  un argument  de sortie pour le s-prog

  USE dimens_m
  USE paramet_m
  USE comgeom
  USE filtreg_v_m, ONLY: filtreg_v
  use rotatf_m, only: rotatf

  IMPLICIT NONE



  ! ..........    variables  en  arguments     .............

  INTEGER, INTENT (IN) :: klevel
  REAL rotin(iim + 1, jjm, klevel), rotout(ip1jm, klevel)

  ! ..........    variables   locales       ................

  REAL ghy(ip1jm, klevel), ghx(ip1jmp1, klevel)
  ! ........................................................


  CALL filtreg_v(rotin, intensive = .true.)

  CALL nxgrad(klevel, rotin, ghx, ghy)
  CALL rotatf(klevel, ghx, ghy, rotout)

  RETURN
END SUBROUTINE laplacien_rot
