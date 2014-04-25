module fxysinus_m

  IMPLICIT NONE

  private
  public fxysinus

contains

  SUBROUTINE fxysinus(rlatu, yprimu, rlatv, yprimv, rlatu1, yprimu1, rlatu2, &
       yprimu2, rlonu, xprimu, rlonv, xprimv, rlonm025, xprimm025, rlonp025, &
       xprimp025)

    ! From LMDZ4/libf/dyn3d/fxysinus.F, version 1.1.1.1, 2004/05/19 12:53:06
    ! and LMDZ4/libf/grid/fxy_sin.h, v 1.1.1.1, 2004/05/19 12:53:05

    ! Author: P. Le Van

    ! Calcul des longitudes et des latitudes pour une fonction f(x, y)
    ! avec y = Arcsin(j).

    USE dimens_m, only: iim, jjm
    USE nr_util, ONLY: pi

    INTEGER i, j

    REAL, intent(out):: rlatu(:), yprimu(:) ! (jjp1)
    REAL, intent(out):: rlatv(:), yprimv(:) ! (jjm)
    REAL, intent(out):: rlatu1(:) ! (jjm)
    real, intent(out):: yprimu1(:), rlatu2(:), yprimu2(:) ! (jjm)
    REAL, intent(out):: rlonu(:), xprimu(:), rlonv(:), xprimv(:) ! (iip1)
    REAL, intent(out):: rlonm025(:) ! (iip1)
    real, intent(out):: xprimm025(:), rlonp025(:), xprimp025(:) ! (iip1)

    ! Local:
    real fxprim

    !-----------------------------------------------------------------------

    fxprim = 2. * pi / iim

    ! Calcul des latitudes et de y' 

    forall(j = 1: jjm + 1)
       rlatu(j) = fy(real(j))
       yprimu(j) = fyprim(real(j))
    END forall

    forall(j = 1: jjm)
       rlatv(j) = fy(real(j) + 0.5)
       rlatu1(j) = fy(real(j) + 0.25)
       rlatu2(j) = fy(real(j) + 0.75)

       yprimv(j) = fyprim(real(j) + 0.5)
       yprimu1(j) = fyprim(real(j) + 0.25)
       yprimu2(j) = fyprim(real(j) + 0.75)
    END forall

    ! Calcul des longitudes et de x' 
    forall(i = 1: iim + 1)
       rlonv(i) = fx(real(i))
       rlonu(i) = fx(real(i) + 0.5)
       rlonm025(i) = fx(real(i) - 0.25)
       rlonp025(i) = fx(real(i) + 0.25)

       xprimv(i) = fxprim
       xprimu(i) = fxprim
       xprimm025(i) = fxprim
       xprimp025(i) = fxprim
    END forall

  END SUBROUTINE fxysinus

  ! Fonctions à changer éventuellement, selon x(x) et y(y) choisis.
  ! Ici, on a l'application particulière suivante :

  ! x = 2 * pi / iim * X
  ! y = pi / jjm * Y

  pure REAL function fy(rj)
    USE dimens_m, only: jjm
    REAL, intent(in):: rj
    fy = asin(1. + 2. * ((1. - rj) / jjm))
  end function fy

  pure real function fx(ri)
    USE dimens_m, only: iim
    USE nr_util, ONLY: pi
    REAL, intent(in):: ri
    fx = 2.*pi/real(iim)*(ri-0.5*real(iim)-1.)
  end function fx

  pure real function fyprim(rj)
    USE dimens_m, only: jjm
    REAL, intent(in):: rj
    fyprim = 1./sqrt((rj-1.)*(jjm + 1.-rj))
  end function fyprim

end module fxysinus_m
