
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/bernoui.F,v 1.1.1.1 2004/05/19
! 12:53:06 lmdzadmin Exp $

SUBROUTINE bernoui(ngrid, nlay, pphi, pecin, pbern)
  USE dimens_m
  USE paramet_m
  USE conf_gcm_m
  USE filtreg_m, ONLY: filtreg
  IMPLICIT NONE

  ! =======================================================================

  ! Auteur:   P. Le Van
  ! -------

  ! Objet:
  ! ------
  ! calcul de la fonction de Bernouilli aux niveaux s  .....
  ! phi  et  ecin  sont des arguments d'entree pour le s-pg .......
  ! bern       est un  argument de sortie pour le s-pg  ......

  ! fonction de Bernouilli = bern = filtre de( geopotentiel +
  ! energ.cinet.)

  ! =======================================================================

  ! -----------------------------------------------------------------------
  ! Decalrations:
  ! -------------


  ! Arguments:
  ! ----------

  INTEGER nlay, ngrid
  REAL, INTENT (IN) :: pphi(ngrid*nlay), pecin(ngrid*nlay)
  REAL pbern(ngrid*nlay)

  ! Local:
  ! ------

  INTEGER ijl

  ! -----------------------------------------------------------------------
  ! calcul de Bernouilli:
  ! ---------------------

  DO ijl = 1, ngrid*nlay
    pbern(ijl) = pphi(ijl) + pecin(ijl)
  END DO

  ! -----------------------------------------------------------------------
  ! filtre:
  ! -------

  CALL filtreg(pbern, jjp1, llm, 2, 1, .TRUE.)

  ! -----------------------------------------------------------------------
  RETURN
END SUBROUTINE bernoui
