!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/gradiv.F,v 1.1.1.1 2004/05/19 12:53:06 lmdzadmin Exp $
!
      SUBROUTINE gradiv(klevel, xcov, ycov, ld, gdx, gdy )
c
c    Auteur :   P. Le Van
c
c   ***************************************************************
c
c                                ld
c       calcul  de  (grad (div) )   du vect. v ....
c
c     xcov et ycov etant les composant.covariantes de v
c   ****************************************************************
c    xcov , ycov et ld  sont des arguments  d'entree pour le s-prog
c     gdx   et  gdy     sont des arguments de sortie pour le s-prog
c
c
      use dimens_m
      use paramet_m
      use logic
            use inidissip_m
      use filtreg_m, only: filtreg
      IMPLICIT NONE
c

      INTEGER klevel
c
      REAL xcov( ip1jmp1,klevel ), ycov( ip1jm,klevel )
      REAL gdx( ip1jmp1,klevel ),   gdy( ip1jm,klevel )

      REAL div(ip1jmp1,llm)

      INTEGER l,ij,iter
      integer, intent(in):: ld
c
c
c
      CALL SCOPY( ip1jmp1*klevel,xcov,1,gdx,1 )
      CALL SCOPY( ip1jm*klevel,  ycov,1,gdy,1 )
c
      DO 10 iter = 1,ld
c
      CALL  diverg( klevel,  gdx , gdy, div          )
      CALL filtreg( div, jjp1, klevel, 2,1, .true.,2 )
      CALL    grad( klevel,  div, gdx, gdy           )
c
      DO 5  l = 1, klevel
      DO 3 ij = 1, ip1jmp1
      gdx( ij,l ) = - gdx( ij,l ) * cdivu
   3  CONTINUE
      DO 4 ij = 1, ip1jm
      gdy( ij,l ) = - gdy( ij,l ) * cdivu
   4  CONTINUE
   5  CONTINUE
c
  10  CONTINUE
      RETURN
      END
