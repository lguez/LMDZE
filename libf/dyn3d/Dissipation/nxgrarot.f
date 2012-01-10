!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/nxgrarot.F,v 1.1.1.1 2004/05/19 12:53:06 lmdzadmin Exp $
!
      SUBROUTINE nxgrarot (klevel,xcov, ycov, lr, grx, gry, crot )
c   ***********************************************************
c
c    Auteur :  P.Le Van  
c
c                                 lr
c      calcul de  ( nXgrad (rot) )   du vect. v  ....
c
c       xcov et ycov  etant les compos. covariantes de  v
c   ***********************************************************
c     xcov , ycov et lr  sont des arguments  d'entree pour le s-prog
c      grx   et  gry     sont des arguments de sortie pour le s-prog
c
c
      use dimens_m
      use paramet_m
      use logic
      use filtreg_m, only: filtreg
      IMPLICIT NONE
c
c
c
      INTEGER klevel
      REAL, intent(in):: xcov( ip1jmp1,klevel ), ycov( ip1jm,klevel )
      REAL, intent(out)::  grx( ip1jmp1,klevel ),  gry( ip1jm,klevel )
      real, intent(in):: crot
c
      REAL rot(ip1jm,llm)

      INTEGER l,ij,iter
      integer, intent(in):: lr
c
c
c
      grx = xcov
      gry = ycov
c
      DO 10 iter = 1,lr
      CALL  rotat (klevel,grx, gry, rot )
      CALL filtreg( rot, jjm, klevel, 2,1, .false.,2)
      CALL nxgrad (klevel,rot, grx, gry )
c
      DO 5  l = 1, klevel
      DO 2 ij = 1, ip1jm
      gry( ij,l ) = - gry( ij,l ) * crot
   2  CONTINUE
      DO 3 ij = 1, ip1jmp1
      grx( ij,l ) = - grx( ij,l ) * crot
   3  CONTINUE
   5  CONTINUE
c
  10  CONTINUE
      RETURN
      END
