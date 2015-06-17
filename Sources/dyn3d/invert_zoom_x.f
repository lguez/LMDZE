module invert_zoom_x_m

  implicit none

  INTEGER, PARAMETER:: nmax = 30000

contains

  subroutine invert_zoom_x(xf, xtild, G, xlon, xprim, xuv)

    use coefpoly_m, only: coefpoly
    USE dimens_m, ONLY: iim
    use dynetat0_m, only: clon
    use nr_util, only: pi_d, twopi_d
    use numer_rec_95, only: hunt

    DOUBLE PRECISION, intent(in):: Xf(0:), xtild(0:), G(0:) ! (0:nmax)

    real, intent(out):: xlon(:), xprim(:) ! (iim)

    DOUBLE PRECISION, intent(in):: xuv
    ! between - 0.25 and 0.5
    ! 0. si calcul aux points scalaires 
    ! 0.5 si calcul aux points U

    ! Local:
    DOUBLE PRECISION Y, abs_y, a0, a1, a2, a3
    integer i, it, iter
    real, parameter:: my_eps = 1e-6

    real xo1, Xf1
    real xvrai(iim), Gvrai(iim) 
    ! intermediary variables because xlon and xprim are simple precision

    !------------------------------------------------------------------

    it = 0 ! initial guess

    DO i = 1, iim
       Y = - pi_d + (i + xuv - 0.75d0) * twopi_d / iim
       ! - pi <= y < pi
       abs_y = abs(y)

       call hunt(xf, abs_y, it, my_lbound = 0)
       ! {0 <= it <= nmax - 1}

       ! Calcul de xvrai(i) et Gvrai(i)

       CALL coefpoly(Xf(it), Xf(it + 1), G(it), G(it + 1), xtild(it), &
            xtild(it + 1), a0, a1, a2, a3)
       xvrai(i) = xtild(it)
       Xf1 = Xf(it)
       Gvrai(i) = G(it)
       xo1 = xvrai(i)
       iter = 1

       do
          xvrai(i) = xvrai(i) - (Xf1 - abs_y) / Gvrai(i)
          IF (ABS(xvrai(i) - xo1) <= my_eps .or. iter == 300) exit
          xo1 = xvrai(i)
          Xf1 = a0 + xvrai(i) * (a1 + xvrai(i) * (a2 + xvrai(i) * a3))
          Gvrai(i) = a1 + xvrai(i) * (2d0 * a2 + xvrai(i) * 3d0 * a3)
       end DO

       if (ABS(xvrai(i) - xo1) > my_eps) then
          ! iter == 300
          print *, 'Pas de solution.'
          print *, i, abs_y
          STOP 1
       end if

       if (y < 0d0) xvrai(i) = - xvrai(i)
    end DO

    DO i = 1, iim -1
       IF (xvrai(i + 1) < xvrai(i)) THEN
          print *, 'xvrai(', i + 1, ') < xvrai(', i, ')'
          STOP 1
       END IF
    END DO

    xlon = xvrai + clon
    xprim = twopi_d / (iim * Gvrai)

  end subroutine invert_zoom_x

end module invert_zoom_x_m
