!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/flumass.F,v 1.1.1.1 2004/05/19 12:53:06 lmdzadmin Exp $
!
      SUBROUTINE flumass (massebx,masseby, vcont, ucont, pbaru, pbarv )

      use dimens_m
      use paramet_m
      use comgeom
      IMPLICIT NONE

c=======================================================================
c
c   Auteurs:  P. Le Van, F. Hourdin  .
c   -------
c
c   Objet:
c   ------
c
c *********************************************************************
c     .... calcul du flux de masse  aux niveaux s ......
c *********************************************************************
c   massebx,masseby,vcont et ucont sont des argum. d'entree pour le s-pg .
c       pbaru  et pbarv            sont des argum.de sortie pour le s-pg .
c
c=======================================================================



      REAL massebx( ip1jmp1,llm ),masseby( ip1jm,llm ) ,
     * vcont( ip1jm,llm ),ucont( ip1jmp1,llm ),pbaru( ip1jmp1,llm ),
     * pbarv( ip1jm,llm )

      REAL apbarun( iip1 ),apbarus( iip1 )

      REAL sairen,saireun,saires,saireus,ctn,cts,ctn0,cts0
      INTEGER  l,ij,i

      REAL       SSUM


      DO  5 l = 1,llm

      DO  1 ij = iip2,ip1jm
      pbaru( ij,l ) = massebx( ij,l ) * ucont( ij,l )
   1  CONTINUE

      DO 3 ij = 1,ip1jm
      pbarv( ij,l ) = masseby( ij,l ) * vcont( ij,l )
   3  CONTINUE

   5  CONTINUE

c    ................................................................
c     calcul de la composante du flux de masse en x aux poles .......
c    ................................................................
c     par la resolution d'1 systeme de 2 equations .

c     la premiere equat.decrivant le calcul de la divergence en 1 point i
c     du pole,ce calcul etant itere de i=1 a i=im .
c                 c.a.d   ,
c     ( ( 0.5*pbaru(i)-0.5*pbaru(i-1) - pbarv(i))/aire(i)   =
c                                           - somme de ( pbarv(n) )/aire pole

c     l'autre equat.specifiant que la moyenne du flux de masse au pole est =0.
c     c.a.d    somme de pbaru(n)*aire locale(n) = 0.

c     on en revient ainsi a determiner la constante additive commune aux pbaru
c     qui representait pbaru(0,j,l) dans l'equat.du calcul de la diverg.au pt
c     i=1 .
c     i variant de 1 a im
c     n variant de 1 a im

      sairen = SSUM( iim,  aire(   1     ), 1 )
      saireun= SSUM( iim, aireu(   1     ), 1 )
      saires = SSUM( iim,  aire( ip1jm+1 ), 1 )
      saireus= SSUM( iim, aireu( ip1jm+1 ), 1 )

      DO 20 l = 1,llm

      ctn =  SSUM( iim, pbarv(    1     ,l),  1 )/ sairen
      cts =  SSUM( iim, pbarv(ip1jmi1+ 1,l),  1 )/ saires

      pbaru(    1   ,l )=   pbarv(    1     ,l ) - ctn * aire(    1    )
      pbaru( ip1jm+1,l )= - pbarv( ip1jmi1+1,l ) + cts * aire( ip1jm+1 )

      DO 11 i = 2,iim
      pbaru(    i    ,l ) = pbaru(   i - 1   ,l )    +
     *                      pbarv(    i      ,l ) - ctn * aire(   i    )

      pbaru( i+ ip1jm,l ) = pbaru( i+ ip1jm-1,l )    -
     *                      pbarv( i+ ip1jmi1,l ) + cts * aire(i+ ip1jm)
  11  CONTINUE
      DO 12 i = 1,iim
      apbarun(i) = aireu(    i   ) * pbaru(   i    , l)
      apbarus(i) = aireu(i +ip1jm) * pbaru(i +ip1jm, l)
  12  CONTINUE
      ctn0 = -SSUM( iim,apbarun,1 )/saireun
      cts0 = -SSUM( iim,apbarus,1 )/saireus
      DO 14 i = 1,iim
      pbaru(   i    , l) = 2. * ( pbaru(   i    , l) + ctn0 )
      pbaru(i+ ip1jm, l) = 2. * ( pbaru(i +ip1jm, l) + cts0 )
  14  CONTINUE

      pbaru(   iip1 ,l ) = pbaru(    1    ,l )
      pbaru( ip1jmp1,l ) = pbaru( ip1jm +1,l )
  20  CONTINUE

      RETURN
      END
