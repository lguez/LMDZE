module invert_zoom_x_m

  implicit none

  DOUBLE PRECISION, private, save:: abs_y
  private funcd

contains

  subroutine invert_zoom_x(beta, xf, xtild, G, xlon, xprim, xuv)

    ! Libraries:
    use nr_util, only: pi_d, twopi_d
    use numer_rec_95, only: hunt, rtsafe    
    
    use coefpoly_m, only: coefpoly, a1, a2, a3
    use dynetat0_chosen_m, only: clon, grossismx
    USE dimensions, ONLY: iim

    DOUBLE PRECISION, intent(in):: beta, Xf(0:), xtild(0:), G(0:) ! (0:nmax)

    real, intent(out):: xlon(:), xprim(:) ! (iim)

    DOUBLE PRECISION, intent(in):: xuv
    ! between - 0.25 and 0.5
    ! 0. si calcul aux points scalaires 
    ! 0.5 si calcul aux points U

    ! Local:
    
    DOUBLE PRECISION Y
    DOUBLE PRECISION h ! step of the uniform grid
    integer i, it

    DOUBLE PRECISION xvrai(iim), Gvrai(iim) 
    ! intermediary variables because xlon and xprim are single precision

    !------------------------------------------------------------------

    print *, "Call sequence information: invert_zoom_x"
    it = 0 ! initial guess
    h = twopi_d / iim

    DO i = 1, iim
       Y = - pi_d + (i + xuv - 0.75d0) * h
       ! - pi <= y < pi
       abs_y = abs(y)

       ! Distinguish boundaries in order to avoid roundoff error.
       ! funcd should be exactly equal to 0 at xtild(it) or xtild(it +
       ! 1) and could be very small with the wrong sign so rtsafe
       ! would fail.
       if (abs_y == 0d0) then
          xvrai(i) = 0d0
          gvrai(i) = grossismx
       else if (abs_y == pi_d) then
          xvrai(i) = pi_d
          gvrai(i) = 2d0 * beta - grossismx
       else
          call hunt(xf, abs_y, it, my_lbound = 0)
          ! {0 <= it <= nmax - 1}

          ! Calcul de xvrai(i) et Gvrai(i)
          CALL coefpoly(Xf(it), Xf(it + 1), G(it), G(it + 1), xtild(it), &
               xtild(it + 1))
          xvrai(i) = rtsafe(funcd, xtild(it), xtild(it + 1), xacc = 1d-6)
          Gvrai(i) = a1 + xvrai(i) * (2d0 * a2 + xvrai(i) * 3d0 * a3)
       end if

       if (y < 0d0) xvrai(i) = - xvrai(i)
    end DO

    DO i = 1, iim - 1
       IF (xvrai(i + 1) < xvrai(i)) THEN
          print *, 'xvrai(', i + 1, ') < xvrai(', i, ')'
          STOP 1
       END IF
    END DO

    xlon = xvrai + clon
    xprim = h / Gvrai

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
