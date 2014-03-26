
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/coefpoly.F,v 1.1.1.1 2004/05/19
! 12:53:05 lmdzadmin Exp $

SUBROUTINE coefpoly(xf1, xf2, xprim1, xprim2, xtild1, xtild2, a0, a1, a2, a3)
  IMPLICIT NONE

  ! ...  Auteur :   P. Le Van  ...


  ! Calcul des coefficients a0, a1, a2, a3 du polynome de degre 3 qui
  ! satisfait aux 4 equations  suivantes :

  ! a0 + a1*xtild1 + a2*xtild1*xtild1 + a3*xtild1*xtild1*xtild1 = Xf1
  ! a0 + a1*xtild2 + a2*xtild2*xtild2 + a3*xtild2*xtild2*xtild2 = Xf2
  ! a1  +     2.*a2*xtild1 +     3.*a3*xtild1*xtild1 = Xprim1
  ! a1  +     2.*a2*xtild2 +     3.*a3*xtild2*xtild2 = Xprim2

  ! On en revient a resoudre un systeme de 4 equat.a 4 inconnues a0,a1,a2,a3

  DOUBLE PRECISION xf1, xf2, xprim1, xprim2, xtild1, xtild2
  DOUBLE PRECISION a1, a2, a3, a0, xtil1car, xtil2car, derr, x1x2car

  xtil1car = xtild1*xtild1
  xtil2car = xtild2*xtild2

  derr = 2.*(xf2-xf1)/(xtild1-xtild2)

  x1x2car = (xtild1-xtild2)*(xtild1-xtild2)

  a3 = (derr+xprim1+xprim2)/x1x2car
  a2 = (xprim1-xprim2+3.*a3*(xtil2car-xtil1car))/(2.*(xtild1-xtild2))

  a1 = xprim1 - 3.*a3*xtil1car - 2.*a2*xtild1
  a0 = xf1 - a3*xtild1*xtil1car - a2*xtil1car - a1*xtild1

END SUBROUTINE coefpoly
