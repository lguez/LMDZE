!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/dudv2.F,v 1.1.1.1 2004/05/19 12:53:06 lmdzadmin Exp $
!
      SUBROUTINE dudv2 ( teta, pkf, bern, du, dv  )

      use dimens_m
      use paramet_m
      use comvert
      IMPLICIT NONE
c
c=======================================================================
c
c   Auteur:  P. Le Van
c   -------
c
c   Objet:
c   ------
c
c   *****************************************************************
c   ..... calcul du terme de pression (gradient de p/densite )   et
c          du terme de ( -gradient de la fonction de Bernouilli ) ...
c   *****************************************************************
c          Ces termes sont ajoutes a  d(ucov)/dt et a d(vcov)/dt  ..
c
c
c    teta , pkf, bern  sont des arguments d'entree  pour le s-pg  ....
c    du et dv          sont des arguments de sortie pour le s-pg  ....
c
c=======================================================================
c

      REAL teta( ip1jmp1,llm ),pkf( ip1jmp1,llm ) ,bern( ip1jmp1,llm ),
     *         du( ip1jmp1,llm ),  dv( ip1jm,llm )
      INTEGER  l,ij
c
c
      DO 5 l = 1,llm
c
      DO 2  ij  = iip2, ip1jm - 1
       du(ij,l) = du(ij,l) + 0.5* ( teta( ij,l ) + teta( ij+1,l ) ) *
     * ( pkf( ij,l ) - pkf(ij+1,l) )  + bern(ij,l) - bern(ij+1,l)
   2  CONTINUE
c
c
c    .....  correction  pour du(iip1,j,l),  j=2,jjm   ......
c    ...          du(iip1,j,l) = du(1,j,l)                 ...
c
CDIR$ IVDEP
      DO 3 ij = iip1+ iip1, ip1jm, iip1
      du( ij,l ) = du( ij - iim,l )
   3  CONTINUE
c
c
      DO 4 ij  = 1,ip1jm
      dv( ij,l) = dv(ij,l) + 0.5 * ( teta(ij,l) + teta( ij+iip1,l ) ) *
     *                             ( pkf(ij+iip1,l) - pkf(  ij,l  ) )
     *                           +   bern( ij+iip1,l ) - bern( ij  ,l )
   4  CONTINUE
c
   5  CONTINUE
c
      RETURN
      END
