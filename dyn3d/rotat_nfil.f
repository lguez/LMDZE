!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/rotat_nfil.F,v 1.1.1.1 2004/05/19 12:53:05 lmdzadmin Exp $
!
      SUBROUTINE rotat_nfil (klevel, x, y, rot )
c
c    Auteur :   P.Le Van 
c**************************************************************
c.          Calcule le rotationnel  non filtre   ,
c      a tous les niveaux d'1 vecteur de comp. x et y ..
c       x  et  y etant des composantes  covariantes  ...
c********************************************************************
c   klevel, x  et y   sont des arguments d'entree pour le s-prog
c        rot          est  un argument  de sortie pour le s-prog
c
      use dimens_m
      use paramet_m
      use comgeom
      IMPLICIT NONE
c
c
c   .....  variables en arguments  ......
c
      INTEGER, intent(in):: klevel
      REAL rot( ip1jm,klevel )
      REAL x( ip1jmp1,klevel ), y( ip1jm,klevel )
c
c  ...   variables  locales  ...
c
      INTEGER  l, ij
c
c
      DO  10 l = 1,klevel
c
        DO   ij = 1, ip1jm - 1
         rot( ij,l )  =    y( ij+1 , l )  -  y( ij,l )   +
     *                   x(ij +iip1, l )  -  x( ij,l )  
        ENDDO
c
c    .... correction pour rot( iip1,j,l)  ....
c    ....   rot(iip1,j,l)= rot(1,j,l) ...
CDIR$ IVDEP
        DO  ij = iip1, ip1jm, iip1
         rot( ij,l ) = rot( ij -iim,l )
        ENDDO
c
  10  CONTINUE

      RETURN
      END
