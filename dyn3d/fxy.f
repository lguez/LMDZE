module fxy_m

  IMPLICIT NONE

  private
  public fxy

contains

  SUBROUTINE fxy(rlatu, yprimu, rlatv, yprimv, rlatu1, yprimu1, rlatu2, &
       yprimu2, rlonu, xprimu, rlonv, xprimv, rlonm025, xprimm025, rlonp025, &
       xprimp025)

    ! From dyn3d/fxy.F, v 1.1.1.1 2004/05/19 12:53:06
    ! Auteur : P. Le Van
    ! Calcul des longitudes et des latitudes pour une fonction f(x, y)
    ! à tangente sinusoïdale et éventuellement avec zoom.

    USE dimens_m, ONLY : iim, jjm

    REAL, INTENT (OUT) :: rlatu(jjm + 1), yprimu(jjm + 1), rlatv(jjm)
    REAL, INTENT (OUT) :: yprimv(jjm)
    REAL, INTENT (OUT) :: rlatu1(jjm)
    REAL, INTENT (OUT) :: yprimu1(jjm), rlatu2(jjm), yprimu2(jjm)
    REAL, INTENT (OUT) :: rlonu(iim + 1), xprimu(iim + 1), rlonv(iim + 1)
    REAL, INTENT (OUT) :: xprimv(iim + 1)
    REAL, INTENT (OUT) :: rlonm025(iim + 1), xprimm025(iim + 1)
    REAL, INTENT (OUT) :: rlonp025(iim + 1)
    REAL, INTENT (OUT) :: xprimp025(iim + 1)

    ! Variables local to the procedure:

    INTEGER i, j

    !------------------------------------------------------------

    ! Calcul des latitudes et de y'

    DO j = 1, jjm + 1
       rlatu(j) = fy(real(j))
       yprimu(j) = fyprim(real(j))
    END DO

    DO j = 1, jjm
       rlatv(j) = fy(real(j) + 0.5)
       rlatu1(j) = fy(real(j) + 0.25)
       rlatu2(j) = fy(real(j) + 0.75)

       yprimv(j) = fyprim(real(j) + 0.5)
       yprimu1(j) = fyprim(real(j) + 0.25)
       yprimu2(j) = fyprim(real(j) + 0.75)
    END DO

    ! Calcul des longitudes et de x'

    DO i = 1, iim + 1
       rlonv(i) = fx(real(i))
       rlonu(i) = fx(real(i) + 0.5)
       rlonm025(i) = fx(real(i) - 0.25)
       rlonp025(i) = fx(real(i) + 0.25)

       xprimv(i) = fxprim(real(i))
       xprimu(i) = fxprim(real(i) + 0.5)
       xprimm025(i) = fxprim(real(i) - 0.25)
       xprimp025(i) = fxprim(real(i) + 0.25)
    END DO

  END SUBROUTINE fxy

  !******************************************************

  ! From grid/fxy_new.h, v 1.1.1.1 2004/05/19 12:53:05

  REAL FUNCTION ripx(ri)

    ! stretching in x

    USE nr_util, ONLY : pi
    USE dimens_m, ONLY : iim

    REAL, INTENT (IN) :: ri

    ripx = (ri - 1.) * 2 * pi / REAL(iim)

  end function ripx

  !******************************************************

  REAL FUNCTION fx(ri)

    ! stretching in x
    USE nr_util, ONLY : pi
    USE serre, ONLY : alphax, pxo, transx
    REAL, INTENT (IN) :: ri

    fx = ripx(ri) + transx + alphax * SIN(ripx(ri) + transx - pxo) - pi

  end function fx

  !******************************************************

  REAL FUNCTION fxprim(ri)

    ! stretching in x

    USE nr_util, ONLY : pi
    USE serre, ONLY : alphax, pxo, transx
    USE dimens_m, ONLY : iim

    REAL, INTENT (IN) :: ri

    fxprim = 2 * pi / REAL(iim) * (1. + alphax * COS(ripx(ri) + transx - pxo))

  end function fxprim

  !******************************************************

  REAL FUNCTION bigy(rj)

    ! stretching in y
    USE nr_util, ONLY : pi
    USE dimens_m, ONLY : jjm
    REAL, INTENT (IN) :: rj

    bigy = 2 * (REAL(jjm + 1) - rj) * pi / jjm

  end function bigy

  !******************************************************

  REAL FUNCTION fy(rj)

    ! stretching in y
    USE nr_util, ONLY : pi
    USE serre, ONLY : alphay, pyo, transy
    REAL, INTENT (IN) :: rj

    fy = (bigy(rj) + transy + alphay * SIN(bigy(rj) + transy - pyo)) / 2 &
         - pi / 2

  end function fy

  !******************************************************

  REAL FUNCTION fyprim(rj)

    ! stretching in y

    USE nr_util, ONLY : pi
    USE serre, ONLY : alphay, pyo, transy
    USE dimens_m, ONLY : jjm

    REAL, INTENT (IN) :: rj

    fyprim = (pi / jjm) * (1. + alphay * COS(bigy(rj) + transy - pyo))

  end function fyprim

end module fxy_m
