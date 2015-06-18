module invert_zoom_x_m

  implicit none

  INTEGER, PARAMETER:: nmax = 30000
  DOUBLE PRECISION abs_y

  private abs_y, funcd

contains

  subroutine invert_zoom_x(xf, xtild, G, xlon, xprim, xuv)

    use coefpoly_m, only: coefpoly, a1, a2, a3
    USE dimens_m, ONLY: iim
    use dynetat0_m, only: clon
    use nr_util, only: pi_d, twopi_d
    use numer_rec_95, only: hunt, rtsafe

    DOUBLE PRECISION, intent(in):: Xf(0:), xtild(0:), G(0:) ! (0:nmax)

    real, intent(out):: xlon(:), xprim(:) ! (iim)

    DOUBLE PRECISION, intent(in):: xuv
    ! between - 0.25 and 0.5
    ! 0. si calcul aux points scalaires 
    ! 0.5 si calcul aux points U

    ! Local:
    DOUBLE PRECISION Y
    integer i, it

    DOUBLE PRECISION xvrai(iim), Gvrai(iim) 
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
            xtild(it + 1))
       xvrai(i) = rtsafe(funcd, xtild(it), xtild(it + 1), xacc = 1d-6)
       Gvrai(i) = a1 + xvrai(i) * (2d0 * a2 + xvrai(i) * 3d0 * a3)
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

  !**********************************************************************

  SUBROUTINE funcd(x, fval, fderiv)

    use coefpoly_m, only: a0, a1, a2, a3

    DOUBLE PRECISION, INTENT(IN):: x
    DOUBLE PRECISION, INTENT(OUT):: fval, fderiv

    fval = a0 + x * (a1 + x * (a2 + x * a3)) - abs_y
    fderiv = a1 + x * (2d0 * a2 + x * 3d0 * a3)

  END SUBROUTINE funcd

end module invert_zoom_x_m
