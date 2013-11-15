!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/coefpoly.F,v 1.1.1.1 2004/05/19 12:53:05 lmdzadmin Exp $
!
      SUBROUTINE coefpoly ( Xf1, Xf2, Xprim1, Xprim2, xtild1,xtild2 ,
     ,                                          a0,a1,a2,a3         )
      IMPLICIT NONE
c
c   ...  Auteur :   P. Le Van  ...
c
c
c    Calcul des coefficients a0, a1, a2, a3 du polynome de degre 3 qui
c      satisfait aux 4 equations  suivantes :

c    a0 + a1*xtild1 + a2*xtild1*xtild1 + a3*xtild1*xtild1*xtild1 = Xf1
c    a0 + a1*xtild2 + a2*xtild2*xtild2 + a3*xtild2*xtild2*xtild2 = Xf2
c               a1  +     2.*a2*xtild1 +     3.*a3*xtild1*xtild1 = Xprim1
c               a1  +     2.*a2*xtild2 +     3.*a3*xtild2*xtild2 = Xprim2

c  On en revient a resoudre un systeme de 4 equat.a 4 inconnues a0,a1,a2,a3

      DOUBLE PRECISION Xf1, Xf2,Xprim1,Xprim2, xtild1,xtild2, xi 
      DOUBLE PRECISION Xfout, Xprim
      DOUBLE PRECISION a1,a2,a3,a0, xtil1car, xtil2car,derr,x1x2car

      xtil1car = xtild1 * xtild1
      xtil2car = xtild2 * xtild2 

      derr= 2. *(Xf2-Xf1)/( xtild1-xtild2)

      x1x2car = ( xtild1-xtild2)*(xtild1-xtild2)

      a3 = (derr + Xprim1+Xprim2 )/x1x2car
      a2     = ( Xprim1 - Xprim2 + 3.* a3 * ( xtil2car-xtil1car ) )    /
     /           (  2.* ( xtild1 - xtild2 )  )

      a1     = Xprim1 -3.* a3 * xtil1car     -2.* a2 * xtild1
      a0     =  Xf1 - a3 * xtild1* xtil1car -a2 * xtil1car - a1 *xtild1

      RETURN
      END
