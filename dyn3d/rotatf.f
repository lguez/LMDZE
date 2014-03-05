
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/rotatf.F,v 1.1.1.1 2004/05/19
! 12:53:05 lmdzadmin Exp $

SUBROUTINE rotatf(klevel, x, y, rot)

  ! Auteur : P.Le Van
  ! **************************************************************
  ! .  calcule le rotationnel
  ! a tous les niveaux d'1 vecteur de comp. x et y ..
  ! x  et  y etant des composantes  covariantes  ...
  ! ********************************************************************
  ! klevel, x  et y   sont des arguments d'entree pour le s-prog
  ! rot          est  un argument  de sortie pour le s-prog

  USE dimens_m
  USE paramet_m
  USE comgeom
  USE filtreg_m, ONLY: filtreg
  IMPLICIT NONE


  ! .....  variables en arguments  ......

  INTEGER, INTENT (IN) :: klevel
  REAL rot(ip1jm, klevel)
  REAL, INTENT (IN) :: x(ip1jmp1, klevel), y(ip1jm, klevel)

  ! ...   variables  locales  ...

  INTEGER l, ij


  DO l = 1, klevel

    DO ij = 1, ip1jm - 1
      rot(ij, l) = y(ij+1, l) - y(ij, l) + x(ij+iip1, l) - x(ij, l)
    END DO

    ! .... correction pour rot( iip1,j,l)  ....
    ! ....   rot(iip1,j,l)= rot(1,j,l) ...
    ! DIR$ IVDEP
    DO ij = iip1, ip1jm, iip1
      rot(ij, l) = rot(ij-iim, l)
    END DO

  END DO

  CALL filtreg(rot, jjm, klevel, 2, 2, .FALSE.)

  DO l = 1, klevel
    DO ij = 1, ip1jm
      rot(ij, l) = rot(ij, l)*unsairez(ij)
    END DO
  END DO


  RETURN
END SUBROUTINE rotatf
