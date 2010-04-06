SUBROUTINE geopot(ngrid,teta,pk,pks,phis,phi)

  ! From libf/dyn3d/geopot.F, version 1.1.1.1 2004/05/19
  ! Author: P. Le Van
  ! Objet : calcul du géopotentiel aux milieux des couches
  ! L'integration se fait de bas en haut.

  USE dimens_m
  USE paramet_m

  IMPLICIT NONE

  ! Arguments:
  INTEGER, INTENT (IN):: ngrid
  REAL, INTENT (IN):: teta(ngrid,llm), pks(ngrid)
  REAL, INTENT (IN) :: phis(ngrid)
  REAL, INTENT (IN) :: pk(ngrid,llm)
  REAL, INTENT (out)::  phi(ngrid,llm)

  ! Local:
  INTEGER l, ij

  ! -----------------------------------------------------------------------

  ! calcul de phi au niveau 1 pres du sol
  DO  ij = 1, ngrid
     phi(ij,1) = phis(ij) + teta(ij,1)*(pks(ij)-pk(ij,1))
  end DO

  ! calcul de phi aux niveaux superieurs
  DO l = 2, llm
     DO ij = 1, ngrid
        phi(ij,l) = phi(ij,l-1) + 0.5 * (teta(ij,l) + teta(ij,l-1)) &
             * (pk(ij,l-1) - pk(ij,l))
     END DO
  END DO

END SUBROUTINE geopot
