!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/nxgraro2.F,v 1.1.1.1 2004/05/19 12:53:06 lmdzadmin Exp $
!
       SUBROUTINE nxgraro2 (klevel,xcov, ycov, lr, grx, gry )
c
c      P.Le Van .
c   ***********************************************************
c                                 lr
c      calcul de  ( nxgrad (rot) )   du vect. v  ....
c
c       xcov et ycov  etant les compos. covariantes de  v
c   ***********************************************************
c     xcov , ycov et lr  sont des arguments  d'entree pour le s-prog
c      grx   et  gry     sont des arguments de sortie pour le s-prog
c
c
      use dimens_m
      use paramet_m
            use inidissip_m
      IMPLICIT NONE
c
c
c    ......  variables en arguments  .......
c
      INTEGER klevel
      REAL xcov( ip1jmp1,klevel ), ycov( ip1jm,klevel )
      REAL  grx( ip1jmp1,klevel ),  gry( ip1jm,klevel )
c
c    ......   variables locales     ........
c
      REAL rot(ip1jm,llm) , signe, nugradrs
      INTEGER l,ij,iter
      integer, intent(in):: lr
c    ........................................................
c
c
c
      signe    = (-1.)**lr
      nugradrs = signe * crot
c
      CALL SCOPY ( ip1jmp1* klevel, xcov, 1, grx, 1 )
      CALL SCOPY (  ip1jm * klevel, ycov, 1, gry, 1 )
c
      CALL     rotatf     ( klevel, grx, gry, rot )
c
      CALL laplacien_rot ( klevel, rot, rot,grx,gry      )

c
c    .....   Iteration de l'operateur laplacien_rotgam  .....
c
      DO  iter = 1, lr -2
        CALL laplacien_rotgam ( klevel, rot, rot )
      ENDDO
c
c
      CALL filtreg( rot, jjm, klevel, 2,1, .FALSE.,1)
      CALL nxgrad ( klevel, rot, grx, gry )
c
      DO    l = 1, klevel
         DO  ij = 1, ip1jm
          gry( ij,l ) = gry( ij,l ) * nugradrs
         ENDDO
         DO  ij = 1, ip1jmp1
          grx( ij,l ) = grx( ij,l ) * nugradrs
         ENDDO
      ENDDO
c
      RETURN
      END
