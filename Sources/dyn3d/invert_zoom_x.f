module invert_zoom_x_m

  implicit none

  INTEGER, PARAMETER:: nmax = 30000

contains

  subroutine invert_zoom_x(xf, xtild, Xprimt, xlon, xprimm, xuv)

    use coefpoly_m, only: coefpoly
    USE dimens_m, ONLY: iim
    use dynetat0_m, only: clon
    use nr_util, only: pi_d, twopi_d
    use numer_rec_95, only: hunt

    DOUBLE PRECISION, intent(in):: Xf(0:), xtild(0:), Xprimt(0:) ! (0:2 * nmax)
    real, intent(out):: xlon(:), xprimm(:) ! (iim)

    DOUBLE PRECISION, intent(in):: xuv
    ! between - 0.25 and 0.5
    ! 0. si calcul aux points scalaires 
    ! 0.5 si calcul aux points U

    ! Local:
    DOUBLE PRECISION xo1, Xfi, a0, a1, a2, a3, Xf1, Xprimin
    integer i, it, iter
    DOUBLE PRECISION, parameter:: my_eps = 1d-6

    DOUBLE PRECISION xxprim(iim), xvrai(iim) 
    ! intermediary variables because xlon and xprimm are simple precision

    !------------------------------------------------------------------

    it = 0 ! initial guess

    DO i = 1, iim
       Xfi = - pi_d + (i + xuv - 0.75d0) * twopi_d / iim
       ! - pi <= xfi < pi

       call hunt(xf, xfi, it, my_lbound = 0)
       it = max(0, it)
       ! In principle, xfi >= xf(0), max with 0 just in case of
       ! roundoff error. {0 <= it <= 2 * nmax - 1}

       ! Calcul de Xf(xvrai(i)) 

       CALL coefpoly(Xf(it), Xf(it + 1), Xprimt(it), Xprimt(it + 1), &
            xtild(it), xtild(it + 1), a0, a1, a2, a3)
       xvrai(i) = xtild(it)
       Xf1 = Xf(it)
       Xprimin = a1 + xvrai(i) * (2d0 * a2 + xvrai(i) * 3d0 * a3)
       xo1 = xvrai(i)
       iter = 1

       do
          xvrai(i) = xvrai(i) - (Xf1 - Xfi) / Xprimin
          IF (ABS(xvrai(i) - xo1) <= my_eps .or. iter == 300) exit
          xo1 = xvrai(i)
          Xf1 = a0 + xvrai(i) * (a1 + xvrai(i) * (a2 + xvrai(i) * a3))
          Xprimin = a1 + xvrai(i) * (2d0 * a2 + xvrai(i) * 3d0 * a3)
       end DO

       if (ABS(xvrai(i) - xo1) > my_eps) then
          ! iter == 300
          print *, 'Pas de solution.'
          print *, i, xfi
          STOP 1
       end if

       xxprim(i) = twopi_d / (iim * Xprimin)
    end DO

    DO i = 1, iim -1
       IF (xvrai(i + 1) < xvrai(i)) THEN
          print *, 'xvrai(', i + 1, ') < xvrai(', i, ')'
          STOP 1
       END IF
    END DO

    xlon = xvrai + clon
    xprimm = xxprim

  end subroutine invert_zoom_x

end module invert_zoom_x_m
