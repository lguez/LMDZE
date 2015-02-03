module fxhyp_loop_ik_m

  implicit none

  INTEGER, PARAMETER:: nmax = 30000

contains

  subroutine fxhyp_loop_ik(ik, decalx, xf, xtild, Xprimt, xzoom, xlon, &
       xprimm, xuv)

    use coefpoly_m, only: coefpoly
    USE dimens_m, ONLY: iim
    use nr_util, only: pi_d, twopi_d, twopi
    use serre, only: grossismx

    INTEGER, intent(in):: ik
    DOUBLE PRECISION, intent(in):: decalx
    DOUBLE PRECISION, intent(in):: Xf(0:), xtild(0:), Xprimt(0:) ! (0:2 * nmax)
    DOUBLE PRECISION, intent(in):: xzoom
    real, intent(out):: xlon(:), xprimm(:) ! (iim + 1)

    DOUBLE PRECISION, intent(in):: xuv
    ! 0. si calcul aux points scalaires 
    ! 0.5 si calcul aux points U 

    ! Local:

    DOUBLE PRECISION xo1, Xfi, xi, a0, a1, a2, a3, Xf1, Xprimin, xi2
    integer ii1, ii2, i, it, iter, idif, ii
    DOUBLE PRECISION, parameter:: my_eps = 1e-3
    DOUBLE PRECISION xxprim(iim + 1), xvrai(iim + 1)
    INTEGER:: is2 = 0

    !------------------------------------------------------------------

    xo1 = 0.

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

       ! Appel de la routine qui calcule les coefficients a0, a1,
       ! a2, a3 d'un polynome de degre 3 qui passe par les points
       ! (Xf(it), xtild(it)) et (Xf(it + 1), xtild(it + 1))

       CALL coefpoly(Xf(it), Xf(it + 1), Xprimt(it), Xprimt(it + 1), &
            xtild(it), xtild(it + 1), a0, a1, a2, a3)

       Xf1 = Xf(it)
       Xprimin = a1 + 2. * a2 * xi + 3. * a3 * xi**2

       iter = 1

       do
          xi = xi - (Xf1 - Xfi) / Xprimin
          IF (ABS(xi - xo1) <= my_eps .or. iter == 300) exit
          xo1 = xi
          xi2 = xi * xi
          Xf1 = a0 + a1 * xi + a2 * xi2 + a3 * xi2 * xi
          Xprimin = a1 + 2. * a2 * xi + 3. * a3 * xi2
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

    DO i = 1, iim
       xlon(i) = xvrai(i)
       xprimm(i) = xxprim(i)
    END DO

    DO i = 1, iim -1
       IF (xvrai(i + 1) < xvrai(i)) THEN
          print *, 'rlonu(', i + 1, ') < rlonu(', i, ')'
          STOP 1
       END IF
    END DO

    IF (ik == 1 .and. (MINval(xvrai(:iim)) < - pi_d - 0.1d0 &
         .or. MAXval(xvrai(:iim)) > pi_d + 0.1d0)) THEN
       IF (xzoom <= 0.) THEN
          i = 1

          do while (xvrai(i) < - pi_d .and. i < iim)
             i = i + 1
          end do

          if (xvrai(i) < - pi_d) then
             print *, 'Xvrai plus petit que - pi !'
             STOP 1
          end if

          is2 = i
       ELSE
          i = iim

          do while (xvrai(i) > pi_d .and. i > 1)
             i = i - 1
          end do

          if (xvrai(i) > pi_d) then
             print *, 'Xvrai plus grand que pi !'
             STOP 1
          end if

          is2 = i
       END IF
    END IF

    if (is2 /= 0) then
       IF (xzoom <= 0.) THEN
          IF (is2 /= 1) THEN
             DO ii = is2, iim
                xlon(ii-is2 + 1) = xvrai(ii)
                xprimm(ii-is2 + 1) = xxprim(ii)
             END DO
             DO ii = 1, is2 -1
                xlon(ii + iim-is2 + 1) = xvrai(ii) + twopi_d
                xprimm(ii + iim-is2 + 1) = xxprim(ii) 
             END DO
          END IF
       else
          idif = iim -is2

          DO ii = 1, is2
             xlon(ii + idif) = xvrai(ii)
             xprimm(ii + idif) = xxprim(ii)
          END DO

          DO ii = 1, idif
             xlon(ii) = xvrai(ii + is2) - twopi_d
             xprimm(ii) = xxprim(ii + is2) 
          END DO
       end IF
    end if

    xlon(iim + 1) = xlon(1) + twopi
    xprimm(iim + 1) = xprimm(1)

  end subroutine fxhyp_loop_ik

end module fxhyp_loop_ik_m
