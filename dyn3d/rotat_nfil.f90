
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/rotat_nfil.F,v 1.1.1.1 2004/05/19
! 12:53:05 lmdzadmin Exp $

SUBROUTINE rotat_nfil(klevel, x, y, rot)

  ! Auteur :   P.Le Van
  ! **************************************************************
  ! .          Calcule le rotationnel  non filtre   ,
  ! a tous les niveaux d'1 vecteur de comp. x et y ..
  ! x  et  y etant des composantes  covariantes  ...
  ! ********************************************************************
  ! klevel, x  et y   sont des arguments d'entree pour le s-prog
  ! rot          est  un argument  de sortie pour le s-prog

  USE dimensions
  USE paramet_m
  USE comgeom
  IMPLICIT NONE


  ! .....  variables en arguments  ......

  INTEGER, INTENT (IN) :: klevel
  REAL rot(ip1jm, klevel)
  REAL x(ip1jmp1, klevel), y(ip1jm, klevel)

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

  RETURN
END SUBROUTINE rotat_nfil
