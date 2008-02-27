!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/vitvert.F,v 1.1.1.1 2004/05/19 12:53:05 lmdzadmin Exp $
!
      SUBROUTINE vitvert ( convm , w )
c

c=======================================================================
c
c   Auteurs:  P. Le Van , F. Hourdin .
c   -------
c
c   Objet:
c   ------
c
c    *******************************************************************
c  .... calcul de la vitesse verticale aux niveaux sigma  ....
c    *******************************************************************
c     convm   est un argument  d'entree pour le s-pg  ......
c       w     est un argument de sortie pour le s-pg  ......
c
c    la vitesse verticale est orientee de  haut en bas .
c    au sol, au niveau sigma(1),   w(i,j,1) = 0.
c    au sommet, au niveau sigma(llm+1) , la vit.verticale est aussi
c    egale a 0. et n'est pas stockee dans le tableau w  .
c
c
c=======================================================================

      use dimens_m
      use paramet_m
      use comvert
      IMPLICIT NONE

      REAL w(ip1jmp1,llm),convm(ip1jmp1,llm)
      INTEGER   l, ij



      DO 2  l = 1,llmm1

      DO 1 ij = 1,ip1jmp1
      w( ij, l+1 ) = convm( ij, l+1 ) - bp(l+1) * convm( ij, 1 )
   1  CONTINUE

   2  CONTINUE

      DO 5 ij  = 1,ip1jmp1
      w(ij,1)  = 0.
5     CONTINUE

      RETURN
      END
