
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/nxgrad_gam.F,v 1.1.1.1 2004/05/19
! 12:53:06 lmdzadmin Exp $

SUBROUTINE nxgrad_gam(klevel, rot, x, y)

  ! P. Le Van

  ! ********************************************************************
  ! calcul du gradient tourne de pi/2 du rotationnel du vect.v
  ! ********************************************************************
  ! rot          est un argument  d'entree pour le s-prog
  ! x  et y    sont des arguments de sortie pour le s-prog

  USE dimens_m
  USE paramet_m
  USE comgeom
  IMPLICIT NONE

  INTEGER, INTENT (IN) :: klevel
  REAL rot(ip1jm, klevel), x(ip1jmp1, klevel), y(ip1jm, klevel)
  INTEGER l, ij

  DO l = 1, klevel

    DO ij = 2, ip1jm
      y(ij, l) = (rot(ij,l)-rot(ij-1,l))*cvscuvgam(ij)
    END DO

    ! ..... correction pour  y ( 1,j,l )  ......

    ! ....    y(1,j,l)= y(iip1,j,l) ....
    ! DIR$ IVDEP
    DO ij = 1, ip1jm, iip1
      y(ij, l) = y(ij+iim, l)
    END DO

    DO ij = iip2, ip1jm
      x(ij, l) = (rot(ij,l)-rot(ij-iip1,l))*cuscvugam(ij)
    END DO
    DO ij = 1, iip1
      x(ij, l) = 0.
      x(ij+ip1jm, l) = 0.
    END DO

  END DO
  RETURN
END SUBROUTINE nxgrad_gam
