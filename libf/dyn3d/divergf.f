!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/divergf.F,v 1.1.1.1 2004/05/19 12:53:05 lmdzadmin Exp $
!
      SUBROUTINE divergf(klevel,x,y,div)
c
c     P. Le Van
c
c  *********************************************************************
c  ... calcule la divergence a tous les niveaux d'1 vecteur de compos. 
c     x et y...
c              x et y  etant des composantes covariantes   ...
c  *********************************************************************
      use dimens_m
      use paramet_m
      use comgeom
      use filtreg_m, only: filtreg
      IMPLICIT NONE
c
c      x  et  y  sont des arguments  d'entree pour le s-prog
c        div      est  un argument  de sortie pour le s-prog
c
c
c   ---------------------------------------------------------------------
c
c    ATTENTION : pendant ce s-pg , ne pas toucher au COMMON/scratch/  .
c
c   ---------------------------------------------------------------------
c
c    ..........          variables en arguments    ...................
c
      INTEGER, intent(in):: klevel
      REAL x( ip1jmp1,klevel ),y( ip1jm,klevel ),div( ip1jmp1,klevel )
      INTEGER   l,ij
c
c    ...............     variables  locales   .........................

      REAL aiy1( iip1 ) , aiy2( iip1 )
      REAL sumypn,sumyps
c    ...................................................................
c
      REAL      SSUM
c
c
      DO 10 l = 1,klevel
c
        DO  ij = iip2, ip1jm - 1
         div( ij + 1, l )     =  
     *   cvusurcu( ij+1 ) * x( ij+1,l ) - cvusurcu( ij ) * x( ij , l) +
     *   cuvsurcv(ij-iim) * y(ij-iim,l) - cuvsurcv(ij+1) * y(ij+1,l) 
        ENDDO
c
c     ....  correction pour  div( 1,j,l)  ......
c     ....   div(1,j,l)= div(iip1,j,l) ....
c
CDIR$ IVDEP
        DO  ij = iip2,ip1jm,iip1
         div( ij,l ) = div( ij + iim,l )
        ENDDO
c
c     ....  calcul  aux poles  .....
c
        DO  ij  = 1,iim
         aiy1(ij) =    cuvsurcv(    ij       ) * y(     ij     , l )
         aiy2(ij) =    cuvsurcv( ij+ ip1jmi1 ) * y( ij+ ip1jmi1, l )
        ENDDO
        sumypn = SSUM ( iim,aiy1,1 ) / apoln
        sumyps = SSUM ( iim,aiy2,1 ) / apols
c
        DO  ij = 1,iip1
         div(     ij    , l ) = - sumypn
         div( ij + ip1jm, l ) =   sumyps
        ENDDO
  10  CONTINUE
c

        CALL filtreg( div, jjp1, klevel, 2, 2, .TRUE., 1 )
      
c
        DO l = 1, klevel
           DO ij = iip2,ip1jm
            div(ij,l) = div(ij,l) * unsaire(ij) 
          ENDDO
        ENDDO
c
       RETURN
       END
