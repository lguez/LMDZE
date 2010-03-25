!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/dteta1.F,v 1.1.1.1 2004/05/19 12:53:06 lmdzadmin Exp $
!
      SUBROUTINE dteta1 ( teta, pbaru, pbarv, dteta)
      use dimens_m
      use paramet_m
      use logic
      use filtreg_m, only: filtreg

      IMPLICIT NONE

c=======================================================================
c
c   Auteur:  P. Le Van
c   -------
c Modif F.Forget 03/94 (on retire q et dq  pour construire dteta1)
c
c   ********************************************************************
c   ... calcul du terme de convergence horizontale du flux d'enthalpie
c        potentielle   ......
c   ********************************************************************
c  .. teta,pbaru et pbarv sont des arguments d'entree  pour le s-pg ....
c     dteta 	          sont des arguments de sortie pour le s-pg ....
c
c=======================================================================



      REAL teta( ip1jmp1,llm ),pbaru( ip1jmp1,llm ),pbarv( ip1jm,llm)
      REAL dteta( ip1jmp1,llm )
      INTEGER   l,ij

      REAL hbyv( ip1jm,llm ), hbxu( ip1jmp1,llm )

c

      DO 5 l = 1,llm

      DO 1  ij = iip2, ip1jm - 1
      hbxu(ij,l) = pbaru(ij,l) * 0.5 * ( teta(ij,l) + teta(ij+1,l) )
   1  CONTINUE

c    .... correction pour  hbxu(iip1,j,l)  .....
c    ....   hbxu(iip1,j,l)= hbxu(1,j,l) ....

CDIR$ IVDEP
      DO 2 ij = iip1+ iip1, ip1jm, iip1
      hbxu( ij, l ) = hbxu( ij - iim, l )
   2  CONTINUE


      DO 3 ij = 1,ip1jm
      hbyv(ij,l)= pbarv(ij, l)* 0.5 * ( teta(ij, l)+ teta(ij +iip1,l) )
   3  CONTINUE

   5  CONTINUE


        CALL  convflu ( hbxu, hbyv, llm, dteta )


c    stockage dans  dh de la convergence horizont. filtree' du  flux
c                  ....                           ...........
c           d'enthalpie potentielle .

      CALL filtreg( dteta, jjp1, llm, 2, 2, .true., 1)

c
      RETURN
      END
