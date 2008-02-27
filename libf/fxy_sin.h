!
! $Header: /home/cvsroot/LMDZ4/libf/grid/fxy_sin.h,v 1.1.1.1 2004/05/19 12:53:05 lmdzadmin Exp $
!
c-----------------------------------------------------------------------
c INCLUDE 'fxyprim.h'
c
c    ................................................................
c    ................  Fonctions in line  ...........................
c    ................................................................
c
      REAL  fy, fx, fxprim, fyprim
      REAL  ri, rj
c
c
      fy(rj)=ASIN(1.+2.*((1.-rj)/FLOAT(jjm)))
      fyprim(rj)=1./SQRT((rj-1.)*(jjm+1.-rj))

      fx    ( ri ) = 2.*pi/FLOAT(iim) * ( ri - 0.5*  FLOAT(iim) - 1. )
c     fx    ( ri ) = 2.*pi/FLOAT(iim) * ( ri - 0.5* ( FLOAT(iim) + 1.) )
      fxprim( ri ) = 2.*pi/FLOAT(iim)
c
c
c    La valeur de pi est passee par le common/const/ou /const2/ .
c    Sinon, il faut la calculer avant d'appeler ces fonctions .
c
c   ----------------------------------------------------------------
c     Fonctions a changer eventuellement, selon x(x) et y(y) choisis .
c   -----------------------------------------------------------------
c
c    .....  ici, on a l'application particuliere suivante   ........
c
c                **************************************
c                **     x = 2. * pi/iim *  X         **
c                **     y =      pi/jjm *  Y         **
c                **************************************
c
c   ..................................................................
c   ..................................................................
c
c
c
c-----------------------------------------------------------------------
