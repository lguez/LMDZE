!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/nxgrad.F,v 1.1.1.1 2004/05/19 12:53:05 lmdzadmin Exp $
!
      SUBROUTINE nxgrad (klevel, rot, x, y )
c
c     P. Le Van
c
c   ********************************************************************
c      calcul du gradient tourne de pi/2 du rotationnel du vect.v
c   ********************************************************************
c       rot          est un argument  d'entree pour le s-prog
c       x  et y    sont des arguments de sortie pour le s-prog
c
      use dimens_m
      use paramet_m
      use comgeom
      IMPLICIT NONE
c
      INTEGER klevel
      REAL rot( ip1jm,klevel ),x( ip1jmp1,klevel ),y(ip1jm,klevel )
      INTEGER   l,ij
c
c
      DO 10 l = 1,klevel
c
      DO 1  ij = 2, ip1jm
      y( ij,l ) = (  rot( ij,l ) - rot( ij-1,l )  ) * cvsurcuv( ij )
   1  CONTINUE
c
c    ..... correction pour  y ( 1,j,l )  ......
c
c    ....    y(1,j,l)= y(iip1,j,l) ....
CDIR$ IVDEP
      DO 2  ij = 1, ip1jm, iip1
      y( ij,l ) = y( ij +iim,l )
   2  CONTINUE
c
      DO 4  ij = iip2,ip1jm
      x( ij,l ) = (  rot( ij,l ) - rot( ij -iip1,l )  ) * cusurcvu( ij )
   4  CONTINUE
      DO 6 ij = 1,iip1
      x(    ij    ,l ) = 0.
      x( ij +ip1jm,l ) = 0.
   6  CONTINUE
c
  10  CONTINUE
      RETURN
      END
