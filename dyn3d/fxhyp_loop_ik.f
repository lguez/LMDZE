module fxhyp_loop_ik_m

  implicit none

  INTEGER, PARAMETER:: nmax = 30000

contains

  subroutine fxhyp_loop_ik(ik, decalx, xf, xtild, Xprimt, xzoom, xlon, &
       xprimm, xuv)

    use coefpoly_m, only: coefpoly
    USE dimens_m, ONLY: iim
    use nr_util, only: pi_d, twopi_d
    use serre, only: grossismx

    INTEGER, intent(in):: ik
    DOUBLE PRECISION, intent(in):: decalx
    DOUBLE PRECISION, intent(in):: Xf(0:), xtild(0:), Xprimt(0:) ! (0:2 * nmax)
    DOUBLE PRECISION, intent(in):: xzoom
    real, intent(out):: xlon(:), xprimm(:) ! (iim)

    DOUBLE PRECISION, intent(in):: xuv
    ! 0. si calcul aux points scalaires 
    ! 0.5 si calcul aux points U 

    ! Local:
    DOUBLE PRECISION xo1, Xfi, xi, a0, a1, a2, a3, Xf1, Xprimin
    integer ii1, ii2, i, it, iter
    DOUBLE PRECISION, parameter:: my_eps = 1e-6
    DOUBLE PRECISION xxprim(iim + 1), xvrai(iim + 1)

    !------------------------------------------------------------------

    IF (ik == 1 .and. grossismx == 1.) THEN
       ii1 = 2 
       ii2 = iim + 1
    else
       ii1=1
       ii2=iim
    END IF

    DO i = ii1, ii2
       Xfi = - pi_d + (i + xuv - decalx) * twopi_d / iim

       it = 2 * nmax
       do while (xfi < xf(it) .and. it >= 1)
          it = it - 1
       end do

       ! Calcul de Xf(xi) 

       xi = xtild(it)

       IF (it == 2 * nmax) THEN
          it = 2 * nmax -1
       END IF

       CALL coefpoly(Xf(it), Xf(it + 1), Xprimt(it), Xprimt(it + 1), &
            xtild(it), xtild(it + 1), a0, a1, a2, a3)
       Xf1 = Xf(it)
       Xprimin = a1 + xi * (2d0 * a2 + xi * 3d0 * a3)
       xo1 = xi
       iter = 1

       do
          xi = xi - (Xf1 - Xfi) / Xprimin
          IF (ABS(xi - xo1) <= my_eps .or. iter == 300) exit
          xo1 = xi
          Xf1 = a0 + xi * (a1 + xi * (a2 + xi * a3))
          Xprimin = a1 + xi * (2d0 * a2 + xi * 3d0 * a3)
       end DO

       if (ABS(xi - xo1) > my_eps) then
          ! iter == 300
          print *, 'Pas de solution.'
          print *, i, xfi
          STOP 1
       end if

       xxprim(i) = twopi_d / (REAL(iim) * Xprimin)
       xvrai(i) = xi + xzoom
    end DO

    IF (ik == 1 .and. grossismx == 1.) THEN
       xvrai(1) = xvrai(iim + 1)-twopi_d
       xxprim(1) = xxprim(iim + 1)
    END IF

    DO i = 1, iim -1
       IF (xvrai(i + 1) < xvrai(i)) THEN
          print *, 'xvrai(', i + 1, ') < xvrai(', i, ')'
          STOP 1
       END IF
    END DO

    DO i = 1, iim
       xlon(i) = xvrai(i)
       xprimm(i) = xxprim(i)
    END DO

  end subroutine fxhyp_loop_ik

end module fxhyp_loop_ik_m
