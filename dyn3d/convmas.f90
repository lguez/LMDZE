
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/convmas.F,v 1.1.1.1 2004/05/19
! 12:53:07 lmdzadmin Exp $

SUBROUTINE convmas(pbaru, pbarv, convm)

  USE dimens_m
  USE paramet_m
  USE disvert_m
  USE conf_gcm_m
  USE filtreg_m, ONLY: filtreg
  IMPLICIT NONE

  ! =======================================================================

  ! Auteurs:  P. Le Van , F. Hourdin  .
  ! -------

  ! Objet:
  ! ------

  ! ********************************************************************
  ! .... calcul de la convergence du flux de masse aux niveaux p ...
  ! ********************************************************************


  ! pbaru  et  pbarv  sont des arguments d'entree pour le s-pg  ....
  ! .....  convm      est  un argument de sortie pour le s-pg  ....

  ! le calcul se fait de haut en bas,
  ! la convergence de masse au niveau p(llm+1) est egale a 0. et
  ! n'est pas stockee dans le tableau convm .


  ! =======================================================================

  ! Declarations:
  ! -------------


  REAL, INTENT (IN) :: pbaru(ip1jmp1, llm), pbarv(ip1jm, llm)
  REAL, INTENT (OUT) :: convm(ip1jmp1, llm)
  INTEGER l, ij


  ! -----------------------------------------------------------------------
  ! ....  calcul de - (d(pbaru)/dx + d(pbarv)/dy ) ......

  CALL convflu(pbaru, pbarv, llm, convm)

  ! -----------------------------------------------------------------------
  ! filtrage:
  ! ---------

  CALL filtreg(convm, jjp1, llm, 2, 2, .TRUE.)

  ! integration de la convergence de masse de haut  en bas ......

  DO l = llmm1, 1, -1
    DO ij = 1, ip1jmp1
      convm(ij, l) = convm(ij, l) + convm(ij, l+1)
    END DO
  END DO

  RETURN
END SUBROUTINE convmas
