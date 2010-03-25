!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/gradiv2.F,v 1.1.1.1 2004/05/19 12:53:07 lmdzadmin Exp $
!
      SUBROUTINE gradiv2(klevel, xcov, ycov, ld, gdx, gdy )
c
c     P. Le Van
c
c   **********************************************************
c                                ld
c       calcul  de  (grad (div) )   du vect. v ....
c
c     xcov et ycov etant les composant.covariantes de v
c   **********************************************************
c     xcont , ycont et ld  sont des arguments  d'entree pour le s-prog
c      gdx   et  gdy       sont des arguments de sortie pour le s-prog
c
c
      use dimens_m
      use paramet_m
      use comgeom
            use inidissip_m
      use filtreg_m, only: filtreg
      IMPLICIT NONE
c
c
c     ........    variables en arguments      ........

      INTEGER klevel
      REAL  xcov( ip1jmp1,klevel ), ycov( ip1jm,klevel )
      REAL   gdx( ip1jmp1,klevel ),  gdy( ip1jm,klevel )
c
c     ........       variables locales       .........
c
      REAL div(ip1jmp1,llm)
      REAL signe, nugrads
      INTEGER l,ij,iter
      integer, intent(in):: ld
      
c    ........................................................
c
c
      CALL SCOPY( ip1jmp1 * klevel, xcov, 1, gdx, 1 )
      CALL SCOPY(   ip1jm * klevel, ycov, 1, gdy, 1 )
c
c
      signe   = (-1.)**ld
      nugrads = signe * cdivu
c


      CALL    divergf( klevel, gdx,   gdy , div )

      IF( ld.GT.1 )   THEN

        CALL laplacien ( klevel, div,  div     )

c    ......  Iteration de l'operateur laplacien_gam   .......

        DO iter = 1, ld -2
         CALL laplacien_gam ( klevel,cuvscvgam1,cvuscugam1,unsair_gam1,
     *                       unsapolnga1, unsapolsga1,  div, div       )
        ENDDO

      ENDIF


       CALL filtreg( div   , jjp1, klevel, 2, 1, .TRUE., 1 )
       CALL  grad  ( klevel,  div,   gdx,  gdy             )

c
       DO   l = 1, klevel
         DO  ij = 1, ip1jmp1
          gdx( ij,l ) = gdx( ij,l ) * nugrads
         ENDDO
         DO  ij = 1, ip1jm
          gdy( ij,l ) = gdy( ij,l ) * nugrads
         ENDDO
       ENDDO
c
       RETURN
       END
