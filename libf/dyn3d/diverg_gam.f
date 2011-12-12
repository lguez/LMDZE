!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/diverg_gam.F,v 1.1.1.1 2004/05/19 12:53:05 lmdzadmin Exp $
!
      SUBROUTINE diverg_gam(klevel,cuvscvgam,cvuscugam,unsairegam ,
     *                       unsapolnga,unsapolsga,  x, y,  div )
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
      REAL cuvscvgam(ip1jm),cvuscugam(ip1jmp1),unsairegam(ip1jmp1)
      REAL unsapolnga,unsapolsga
c
c    ...............     variables  locales   .........................

      REAL aiy1( iip1 ) , aiy2( iip1 )
      REAL sumypn,sumyps
      INTEGER   l,ij
c    ...................................................................
c
      REAL      SSUM
c
c
      DO 10 l = 1,klevel
c
        DO  ij = iip2, ip1jm - 1
         div( ij + 1, l )     = (  
     *  cvuscugam( ij+1 ) * x( ij+1,l ) - cvuscugam( ij ) * x( ij , l) +
     *  cuvscvgam(ij-iim) * y(ij-iim,l) - cuvscvgam(ij+1) * y(ij+1,l) )* 
     *         unsairegam( ij+1 )
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
         aiy1(ij) =    cuvscvgam(    ij       ) * y(     ij     , l )
         aiy2(ij) =    cuvscvgam( ij+ ip1jmi1 ) * y( ij+ ip1jmi1, l )
        ENDDO
        sumypn = SSUM ( iim,aiy1,1 ) * unsapolnga
        sumyps = SSUM ( iim,aiy2,1 ) * unsapolsga
c
        DO  ij = 1,iip1
         div(     ij    , l ) = - sumypn 
         div( ij + ip1jm, l ) =   sumyps 
        ENDDO
  10  CONTINUE
c

       RETURN
       END
