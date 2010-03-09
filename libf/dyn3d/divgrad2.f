!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/divgrad2.F,v 1.1.1.1 2004/05/19 12:53:06 lmdzadmin Exp $
!
      SUBROUTINE divgrad2 ( klevel, h, deltapres, lh, divgra )
c
c     P. Le Van
c
c   ***************************************************************
c
c     .....   calcul de  (div( grad ))   de (  pext * h ) .....
c   ****************************************************************
c   h ,klevel,lh et pext  sont des arguments  d'entree pour le s-prg
c         divgra     est  un argument  de sortie pour le s-prg
c
      use dimens_m
      use paramet_m
      use comgeom
            use inidissip_m
      IMPLICIT NONE
c

c    .......    variables en arguments   .......
c
      INTEGER klevel
      REAL h( ip1jmp1,klevel ), deltapres( ip1jmp1,klevel )
      REAL divgra( ip1jmp1,klevel)
c
c    .......    variables  locales    ..........
c
      REAL     signe, nudivgrs, sqrtps( ip1jmp1,llm )
      INTEGER  l,ij,iter
      integer, intent(in):: lh
c    ...................................................................

c
      signe    = (-1.)**lh
      nudivgrs = signe * cdivh

      CALL SCOPY ( ip1jmp1 * klevel, h, 1, divgra, 1 )

c
      CALL laplacien( klevel, divgra, divgra )
     
      DO l = 1, klevel
       DO ij = 1, ip1jmp1
        sqrtps( ij,l ) = SQRT( deltapres(ij,l) )
       ENDDO
      ENDDO
c
      DO l = 1, klevel
        DO ij = 1, ip1jmp1
         divgra(ij,l) = divgra(ij,l) * sqrtps(ij,l)
        ENDDO
      ENDDO
   
c    ........    Iteration de l'operateur  laplacien_gam    ........
c
      DO  iter = 1, lh - 2
       CALL laplacien_gam ( klevel,cuvscvgam2,cvuscugam2,unsair_gam2,
     *                     unsapolnga2, unsapolsga2,  divgra, divgra )
      ENDDO
c
c    ...............................................................
 
      DO l = 1, klevel
        DO ij = 1, ip1jmp1
          divgra(ij,l) = divgra(ij,l) * sqrtps(ij,l)
        ENDDO
      ENDDO
c
      CALL laplacien ( klevel, divgra, divgra )
c
      DO l  = 1,klevel
      DO ij = 1,ip1jmp1
      divgra(ij,l) =  nudivgrs * divgra(ij,l) / deltapres(ij,l)
      ENDDO
      ENDDO

      RETURN
      END
