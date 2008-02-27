!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/convflu.F,v 1.1.1.1 2004/05/19 12:53:05 lmdzadmin Exp $
!
      SUBROUTINE convflu( xflu,yflu,nbniv,convfl )
c
c  P. Le Van
c
c
c    *******************************************************************
c  ... calcule la (convergence horiz. * aire locale)du flux ayant pour
c      composantes xflu et yflu ,variables extensives .  ......
c    *******************************************************************
c      xflu , yflu et nbniv sont des arguments d'entree pour le s-pg ..
c      convfl                est  un argument de sortie pour le s-pg .
c
c     njxflu  est le nombre de lignes de latitude de xflu, 
c     ( = jjm ou jjp1 )
c     nbniv   est le nombre de niveaux vert. de  xflu et de yflu .
c
      use dimens_m
      use paramet_m
      use comgeom
      IMPLICIT NONE
c
      REAL       xflu,yflu,convfl,convpn,convps
      INTEGER    l,ij,nbniv
      DIMENSION  xflu( ip1jmp1,nbniv ),yflu( ip1jm,nbniv ) ,
     *         convfl( ip1jmp1,nbniv )
c
      REAL       SSUM
c
c
c
      DO 5 l = 1,nbniv
c
      DO 2  ij = iip2, ip1jm - 1
      convfl( ij + 1,l ) =  xflu(   ij,l ) - xflu( ij +  1,l )   +
     *                      yflu(ij +1,l ) - yflu( ij -iim,l )
   2  CONTINUE
c
c

c     ....  correction pour  convfl( 1,j,l)  ......
c     ....   convfl(1,j,l)= convfl(iip1,j,l) ...
c
CDIR$ IVDEP
      DO 3 ij = iip2,ip1jm,iip1
      convfl( ij,l ) = convfl( ij + iim,l )
   3  CONTINUE
c
c     ......  calcul aux poles  .......
c
      convpn =   SSUM( iim, yflu(     1    ,l ),  1 )
      convps = - SSUM( iim, yflu( ip1jm-iim,l ),  1 )
      DO 4 ij = 1,iip1
      convfl(     ij   ,l ) = convpn * aire(   ij     ) / apoln
      convfl( ij+ ip1jm,l ) = convps * aire( ij+ ip1jm) / apols
   4  CONTINUE
c
   5  CONTINUE
      RETURN
      END
