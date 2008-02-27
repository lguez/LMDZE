!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/massbarxy.F,v 1.1.1.1 2004/05/19 12:53:07 lmdzadmin Exp $
!
      SUBROUTINE massbarxy(  masse, massebxy )
c
c **********************************************************************
c
c  Calcule les moyennes en x et  y de la masse d'air dans chaque maille.
c **********************************************************************
c    Auteurs : P. Le Van , Fr. Hourdin  .
c   ..........
c
c  ..  masse          est  un  argum. d'entree  pour le s-pg ...
c  ..  massebxy       est  un  argum. de sortie pour le s-pg ...
c     
c
c     IMPLICIT NONE
c
      use dimens_m
      use paramet_m
      use comconst
      use comgeom
c
       REAL  masse( ip1jmp1,llm ), massebxy( ip1jm,llm )
c

      DO   100    l = 1 , llm
c
      DO 5 ij = 1, ip1jm - 1
      massebxy( ij,l ) = masse(    ij  ,l ) * alpha2(   ij    )   +
     +                   masse(   ij+1 ,l ) * alpha3(  ij+1   )   +
     +                   masse( ij+iip1,l ) * alpha1( ij+iip1 )   +
     +                   masse( ij+iip2,l ) * alpha4( ij+iip2 )
   5  CONTINUE

c    ....  correction pour     massebxy( iip1,j )  ........

CDIR$ IVDEP

      DO 7 ij = iip1, ip1jm, iip1
      massebxy( ij,l ) = massebxy( ij - iim,l )
   7  CONTINUE

100   CONTINUE
c
      RETURN
      END
