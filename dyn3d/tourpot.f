!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/tourpot.F,v 1.1.1.1 2004/05/19 12:53:06 lmdzadmin Exp $
!
      SUBROUTINE tourpot ( vcov, ucov, massebxy, vorpot )

c=======================================================================
c
c   Auteur:  P. Le Van
c   -------
c
c   Objet:
c   ------
c
c    *******************************************************************
c    .........      calcul du tourbillon potentiel             .........
c    *******************************************************************
c
c     vcov,ucov,fext et pbarxyfl sont des argum. d'entree pour le s-pg .
c             vorpot            est  un argum.de sortie pour le s-pg .
c
c=======================================================================

      use dimens_m
      use paramet_m
      use conf_gcm_m
      use comgeom
      use filtreg_m, only: filtreg

      IMPLICIT NONE

      REAL  rot( ip1jm,llm )
      REAL, intent(in):: vcov( ip1jm,llm ),ucov( ip1jmp1,llm )
      REAL massebxy( ip1jm,llm ),vorpot( ip1jm,llm )

      INTEGER l, ij




c  ... vorpot = ( Filtre( d(vcov)/dx - d(ucov)/dy ) + fext ) /psbarxy ..



c    ........  Calcul du rotationnel du vent V  puis filtrage  ........

      DO 5 l = 1,llm

      DO 2 ij = 1, ip1jm - 1
      rot( ij,l ) = vcov(ij+1,l)-vcov(ij,l)+ucov(ij+iip1,l)-ucov(ij,l)
   2  CONTINUE

c    ....  correction pour  rot( iip1,j,l )  .....
c    ....     rot(iip1,j,l) = rot(1,j,l)    .....

CDIR$ IVDEP

      DO 3 ij = iip1, ip1jm, iip1
      rot( ij,l ) = rot( ij -iim, l )
   3  CONTINUE

   5  CONTINUE


      CALL  filtreg( rot, jjm, llm, 2, 1, .FALSE.)


      DO 10 l = 1, llm

      DO 6 ij = 1, ip1jm - 1
      vorpot( ij,l ) = ( rot(ij,l) + fext(ij) ) / massebxy(ij,l)
   6  CONTINUE

c    ..... correction pour  vorpot( iip1,j,l)  .....
c    ....   vorpot(iip1,j,l)= vorpot(1,j,l) ....
CDIR$ IVDEP
      DO 8 ij = iip1, ip1jm, iip1
      vorpot( ij,l ) = vorpot( ij -iim,l )
   8  CONTINUE

  10  CONTINUE

      RETURN
      END
