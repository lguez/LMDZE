
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/dudv1.F,v 1.1.1.1 2004/05/19
! 12:53:06 lmdzadmin Exp $

SUBROUTINE dudv1(vorpot, pbaru, pbarv, du, dv)
  USE dimens_m
  USE paramet_m
  IMPLICIT NONE

  ! -----------------------------------------------------------------------

  ! Auteur:   P. Le Van
  ! -------

  ! Objet:
  ! ------
  ! calcul du terme de  rotation
  ! ce terme est ajoute a  d(ucov)/dt et a d(vcov)/dt  ..
  ! vorpot, pbaru et pbarv sont des arguments d'entree  pour le s-pg ..
  ! du  et dv              sont des arguments de sortie pour le s-pg ..

  ! -----------------------------------------------------------------------


  REAL vorpot(ip1jm, llm), pbaru(ip1jmp1, llm), pbarv(ip1jm, llm), &
    du(ip1jmp1, llm), dv(ip1jm, llm)
  INTEGER l, ij


  DO l = 1, llm

    DO ij = iip2, ip1jm - 1
      du(ij, l) = 0.125*(vorpot(ij-iip1,l)+vorpot(ij,l))* &
        (pbarv(ij-iip1,l)+pbarv(ij-iim,l)+pbarv(ij,l)+pbarv(ij+1,l))
    END DO

    DO ij = 1, ip1jm - 1
      dv(ij+1, l) = -0.125*(vorpot(ij,l)+vorpot(ij+1,l))*(pbaru(ij,l)+pbaru( &
        ij+1,l)+pbaru(ij+iip1,l)+pbaru(ij+iip2,l))
    END DO

    ! .... correction  pour  dv( 1,j,l )  .....
    ! ....   dv(1,j,l)= dv(iip1,j,l) ....

    ! DIR$ IVDEP
    DO ij = 1, ip1jm, iip1
      dv(ij, l) = dv(ij+iim, l)
    END DO

  END DO
  RETURN
END SUBROUTINE dudv1
