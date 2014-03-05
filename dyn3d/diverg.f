
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/diverg.F,v 1.1.1.1 2004/05/19
! 12:53:06 lmdzadmin Exp $

SUBROUTINE diverg(klevel, x, y, div)

  ! P. Le Van

  ! *********************************************************************
  ! ... calcule la divergence a tous les niveaux d'1 vecteur de compos.
  ! x et y...
  ! x et y  etant des composantes covariantes   ...
  ! *********************************************************************
  USE dimens_m
  USE paramet_m
  USE comgeom
  IMPLICIT NONE

  ! x  et  y  sont des arguments  d'entree pour le s-prog
  ! div      est  un argument  de sortie pour le s-prog


  ! ---------------------------------------------------------------------

  ! ATTENTION : pendant ce s-pg , ne pas toucher au COMMON/scratch/  .

  ! ---------------------------------------------------------------------

  ! ..........          variables en arguments    ...................

  INTEGER klevel
  REAL x(ip1jmp1, klevel), y(ip1jm, klevel), div(ip1jmp1, klevel)
  INTEGER l, ij

  ! ...............     variables  locales   .........................

  REAL aiy1(iip1), aiy2(iip1)
  REAL sumypn, sumyps
  ! ...................................................................

  REAL ssum


  DO l = 1, klevel

    DO ij = iip2, ip1jm - 1
      div(ij+1, l) = cvusurcu(ij+1)*x(ij+1, l) - cvusurcu(ij)*x(ij, l) + &
        cuvsurcv(ij-iim)*y(ij-iim, l) - cuvsurcv(ij+1)*y(ij+1, l)
    END DO

    ! ....  correction pour  div( 1,j,l)  ......
    ! ....   div(1,j,l)= div(iip1,j,l) ....

    ! DIR$ IVDEP
    DO ij = iip2, ip1jm, iip1
      div(ij, l) = div(ij+iim, l)
    END DO

    ! ....  calcul  aux poles  .....

    DO ij = 1, iim
      aiy1(ij) = cuvsurcv(ij)*y(ij, l)
      aiy2(ij) = cuvsurcv(ij+ip1jmi1)*y(ij+ip1jmi1, l)
    END DO
    sumypn = ssum(iim, aiy1, 1)/apoln
    sumyps = ssum(iim, aiy2, 1)/apols

    DO ij = 1, iip1
      div(ij, l) = -sumypn
      div(ij+ip1jm, l) = sumyps
    END DO
  END DO


  ! cc        CALL filtreg( div, jjp1, klevel, 2, 2, .TRUE., 1 )


  DO l = 1, klevel
    DO ij = iip2, ip1jm
      div(ij, l) = div(ij, l)*unsaire(ij)
    END DO
  END DO

  RETURN
END SUBROUTINE diverg
