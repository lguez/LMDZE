
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/covcont.F,v 1.1.1.1 2004/05/19
! 12:53:07 lmdzadmin Exp $

SUBROUTINE covcont(klevel, ucov, vcov, ucont, vcont)
  USE dimens_m
  USE paramet_m
  USE comgeom
  IMPLICIT NONE

  ! =======================================================================

  ! Auteur:  P. Le Van
  ! -------

  ! Objet:
  ! ------

  ! *********************************************************************
  ! calcul des compos. contravariantes a partir des comp.covariantes
  ! ********************************************************************

  ! =======================================================================


  INTEGER klevel
  REAL, INTENT (IN) :: ucov(ip1jmp1, klevel), vcov(ip1jm, klevel)
  REAL ucont(ip1jmp1, klevel), vcont(ip1jm, klevel)
  INTEGER l, ij


  DO l = 1, klevel

    DO ij = iip2, ip1jm
      ucont(ij, l) = ucov(ij, l)*unscu2(ij)
    END DO

    DO ij = 1, ip1jm
      vcont(ij, l) = vcov(ij, l)*unscv2(ij)
    END DO

  END DO
  RETURN
END SUBROUTINE covcont
