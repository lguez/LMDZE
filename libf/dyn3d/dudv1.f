!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/dudv1.F,v 1.1.1.1 2004/05/19 12:53:06 lmdzadmin Exp $
!
      SUBROUTINE dudv1 ( vorpot, pbaru, pbarv, du, dv )
      use dimens_m
      use paramet_m
      IMPLICIT NONE
c
c-----------------------------------------------------------------------
c
c   Auteur:   P. Le Van
c   -------
c
c   Objet:
c   ------
c   calcul du terme de  rotation
c   ce terme est ajoute a  d(ucov)/dt et a d(vcov)/dt  ..
c   vorpot, pbaru et pbarv sont des arguments d'entree  pour le s-pg ..
c   du  et dv              sont des arguments de sortie pour le s-pg ..
c
c-----------------------------------------------------------------------


      REAL vorpot( ip1jm,llm ) ,pbaru( ip1jmp1,llm ) ,
     *     pbarv( ip1jm,llm ) ,du( ip1jmp1,llm ) ,dv( ip1jm,llm )
      INTEGER  l,ij
c
c
      DO 10 l = 1,llm
c
      DO 2  ij = iip2, ip1jm - 1
      du( ij,l ) = 0.125 *(  vorpot(ij-iip1, l) + vorpot( ij, l)  ) *
     *                    (   pbarv(ij-iip1, l) + pbarv(ij-iim,  l) +
     *                        pbarv(   ij  , l) + pbarv(ij+ 1 ,  l)   )
   2  CONTINUE
c
      DO 3 ij = 1, ip1jm - 1
      dv( ij+1,l ) = - 0.125 *(  vorpot(ij, l)  + vorpot(ij+1, l)  ) *
     *                        (   pbaru(ij, l)  +  pbaru(ij+1   , l) +
     *                       pbaru(ij+iip1, l)  +  pbaru(ij+iip2, l)  )
   3  CONTINUE
c
c    .... correction  pour  dv( 1,j,l )  .....
c    ....   dv(1,j,l)= dv(iip1,j,l) ....
c
CDIR$ IVDEP
      DO 4 ij = 1, ip1jm, iip1
      dv( ij,l ) = dv( ij + iim, l )
   4  CONTINUE
c
  10  CONTINUE
      RETURN
      END
