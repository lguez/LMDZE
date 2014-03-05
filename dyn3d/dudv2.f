
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/dudv2.F,v 1.1.1.1 2004/05/19
! 12:53:06 lmdzadmin Exp $

SUBROUTINE dudv2(teta, pkf, bern, du, dv)

  USE dimens_m
  USE paramet_m
  USE disvert_m
  IMPLICIT NONE

  ! =======================================================================

  ! Auteur:  P. Le Van
  ! -------

  ! Objet:
  ! ------

  ! *****************************************************************
  ! ..... calcul du terme de pression (gradient de p/densite )   et
  ! du terme de ( -gradient de la fonction de Bernouilli ) ...
  ! *****************************************************************
  ! Ces termes sont ajoutes a  d(ucov)/dt et a d(vcov)/dt  ..


  ! teta , pkf, bern  sont des arguments d'entree  pour le s-pg  ....
  ! du et dv          sont des arguments de sortie pour le s-pg  ....

  ! =======================================================================


  REAL, INTENT (IN) :: teta(ip1jmp1, llm)
  REAL pkf(ip1jmp1, llm), bern(ip1jmp1, llm), du(ip1jmp1, llm), &
    dv(ip1jm, llm)
  INTEGER l, ij


  DO l = 1, llm

    DO ij = iip2, ip1jm - 1
      du(ij, l) = du(ij, l) + 0.5*(teta(ij,l)+teta(ij+1,l))*(pkf(ij,l)-pkf(ij &
        +1,l)) + bern(ij, l) - bern(ij+1, l)
    END DO


    ! .....  correction  pour du(iip1,j,l),  j=2,jjm   ......
    ! ...          du(iip1,j,l) = du(1,j,l)                 ...

    ! DIR$ IVDEP
    DO ij = iip1 + iip1, ip1jm, iip1
      du(ij, l) = du(ij-iim, l)
    END DO


    DO ij = 1, ip1jm
      dv(ij, l) = dv(ij, l) + 0.5*(teta(ij,l)+teta(ij+iip1,l))*(pkf(ij+iip1,l &
        )-pkf(ij,l)) + bern(ij+iip1, l) - bern(ij, l)
    END DO

  END DO

  RETURN
END SUBROUTINE dudv2
