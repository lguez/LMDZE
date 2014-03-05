
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/fxysinus.F,v 1.1.1.1 2004/05/19
! 12:53:06 lmdzadmin Exp $

SUBROUTINE fxysinus(rlatu, yprimu, rlatv, yprimv, rlatu1, yprimu1, rlatu2, &
    yprimu2, rlonu, xprimu, rlonv, xprimv, rlonm025, xprimm025, rlonp025, &
    xprimp025)


  USE dimens_m
  USE paramet_m
  USE comconst
  USE nr_util, ONLY: pi
  IMPLICIT NONE

  ! Calcul  des longitudes et des latitudes  pour une fonction f(x,y)
  ! avec y = Asin( j )  .

  ! Auteur  :  P. Le Van



  INTEGER i, j

  REAL rlatu(jjp1), yprimu(jjp1), rlatv(jjm), yprimv(jjm), rlatu1(jjm), &
    yprimu1(jjm), rlatu2(jjm), yprimu2(jjm)
  REAL rlonu(iip1), xprimu(iip1), rlonv(iip1), xprimv(iip1), rlonm025(iip1), &
    xprimm025(iip1), rlonp025(iip1), xprimp025(iip1)


  ! $Header: /home/cvsroot/LMDZ4/libf/grid/fxy_sin.h,v 1.1.1.1 2004/05/19
  ! 12:53:05 lmdzadmin Exp $

  ! -----------------------------------------------------------------------

  ! ................................................................
  ! ................  Fonctions in line  ...........................
  ! ................................................................

  REAL fy, fx, fxprim, fyprim
  REAL ri, rj


  fy(rj) = asin(1.+2.*((1.-rj)/float(jjm)))
  fyprim(rj) = 1./sqrt((rj-1.)*(jjm+1.-rj))

  fx(ri) = 2.*pi/float(iim)*(ri-0.5*float(iim)-1.)
  ! fx    ( ri ) = 2.*pi/FLOAT(iim) * ( ri - 0.5* ( FLOAT(iim) + 1.) )
  fxprim(ri) = 2.*pi/float(iim)


  ! La valeur de pi est passee par le common/const/ou /const2/ .
  ! Sinon, il faut la calculer avant d'appeler ces fonctions .

  ! ----------------------------------------------------------------
  ! Fonctions a changer eventuellement, selon x(x) et y(y) choisis .
  ! -----------------------------------------------------------------

  ! .....  ici, on a l'application particuliere suivante   ........

  ! **************************************
  ! **     x = 2. * pi/iim *  X         **
  ! **     y =      pi/jjm *  Y         **
  ! **************************************

  ! ..................................................................
  ! ..................................................................



  ! -----------------------------------------------------------------------

  ! ......  calcul  des  latitudes  et de y'   .....

  DO j = 1, jjm + 1
    rlatu(j) = fy(float(j))
    yprimu(j) = fyprim(float(j))
  END DO


  DO j = 1, jjm

    rlatv(j) = fy(float(j)+0.5)
    rlatu1(j) = fy(float(j)+0.25)
    rlatu2(j) = fy(float(j)+0.75)

    yprimv(j) = fyprim(float(j)+0.5)
    yprimu1(j) = fyprim(float(j)+0.25)
    yprimu2(j) = fyprim(float(j)+0.75)

  END DO


  ! .....  calcul   des  longitudes et de  x'   .....

  DO i = 1, iim + 1
    rlonv(i) = fx(float(i))
    rlonu(i) = fx(float(i)+0.5)
    rlonm025(i) = fx(float(i)-0.25)
    rlonp025(i) = fx(float(i)+0.25)

    xprimv(i) = fxprim(float(i))
    xprimu(i) = fxprim(float(i)+0.5)
    xprimm025(i) = fxprim(float(i)-0.25)
    xprimp025(i) = fxprim(float(i)+0.25)
  END DO


  RETURN
END SUBROUTINE fxysinus

