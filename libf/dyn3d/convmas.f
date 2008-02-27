!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/convmas.F,v 1.1.1.1 2004/05/19 12:53:07 lmdzadmin Exp $
!
      SUBROUTINE convmas (pbaru, pbarv, convm )
c
      use dimens_m
      use paramet_m
      use comvert
      use logic
      IMPLICIT NONE

c=======================================================================
c
c   Auteurs:  P. Le Van , F. Hourdin  .
c   -------
c
c   Objet:
c   ------
c
c   ********************************************************************
c   .... calcul de la convergence du flux de masse aux niveaux p ...
c   ********************************************************************
c
c
c     pbaru  et  pbarv  sont des arguments d'entree pour le s-pg  ....
c      .....  convm      est  un argument de sortie pour le s-pg  ....
c
c    le calcul se fait de haut en bas, 
c    la convergence de masse au niveau p(llm+1) est egale a 0. et
c    n'est pas stockee dans le tableau convm .
c
c
c=======================================================================
c
c   Declarations:
c   -------------


      REAL pbaru( ip1jmp1,llm ),pbarv( ip1jm,llm ),convm(  ip1jmp1,llm )
      INTEGER   l,ij


c-----------------------------------------------------------------------
c    ....  calcul de - (d(pbaru)/dx + d(pbarv)/dy ) ......

      CALL  convflu( pbaru, pbarv, llm, convm )

c-----------------------------------------------------------------------
c   filtrage:
c   ---------

       CALL filtreg( convm, jjp1, llm, 2, 2, .true., 1 )

c    integration de la convergence de masse de haut  en bas ......

      DO      l      = llmm1, 1, -1
        DO    ij     = 1, ip1jmp1
         convm(ij,l) = convm(ij,l) + convm(ij,l+1)
        ENDDO
      ENDDO
c
      RETURN
      END
