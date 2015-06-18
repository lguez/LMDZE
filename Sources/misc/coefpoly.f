module coefpoly_m

  IMPLICIT NONE

  DOUBLE PRECISION a0, a1, a2, a3

contains

  SUBROUTINE coefpoly(y1, y2, yp1, yp2, x1, x2)

    ! From LMDZ4/libf/dyn3d/coefpoly.F, version 1.1.1.1 2004/05/19 12:53:05

    ! Author: P. Le Van

    ! Calcul des coefficients a0, a1, a2, a3 du polynôme de degré 3
    ! qui passe par les points (x1, Y1) et (x2, Y2) avec les
    ! dérivées yp1 et yp2. Système linéaire de 4 équations à 4
    ! inconnues :

    ! a0 + a1 * x1 + a2 * x1**2 + a3 * x1**3 = Y1
    ! a0 + a1 * x2 + a2 * x2**2 + a3 * x2**3 = Y2
    ! a1 + 2 * a2 * x1 + 3 * a3 * x1**2 = Yp1
    ! a1 + 2 * a2 * x2 + 3 * a3 * x2**2 = Yp2

    DOUBLE PRECISION, intent(in):: y1, y2, yp1, yp2, x1, x2

    ! Local:
    DOUBLE PRECISION x1car, x2car

    !------------------------------------------------------------

    x1car = x1 * x1
    x2car = x2 * x2

    a3 = (2d0 * (y2-y1)/(x1-x2)+yp1+yp2)/((x1-x2) * (x1-x2))
    a2 = (yp1-yp2+3d0 * a3 * (x2car-x1car))/(2d0 * (x1-x2))

    a1 = yp1 - 3d0 * a3 * x1car - 2d0 * a2 * x1
    a0 = y1 - a3 * x1 * x1car - a2 * x1car - a1 * x1

  END SUBROUTINE coefpoly

end module coefpoly_m
