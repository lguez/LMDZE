module fxhyp_m

  IMPLICIT NONE

contains

  SUBROUTINE fxhyp(xprimm025, rlonv, xprimv, rlonu, xprimu, xprimp025)

    ! From LMDZ4/libf/dyn3d/fxhyp.F, version 1.2, 2005/06/03 09:11:32
    ! Author: P. Le Van, from formulas by R. Sadourny

    ! Calcule les longitudes et dérivées dans la grille du GCM pour
    ! une fonction f(x) à dérivée tangente hyperbolique.

    ! Il vaut mieux avoir : grossismx \times dzoom < pi

    ! Le premier point scalaire pour une grille regulière (grossismx =
    ! 1., taux = 0., clon = 0.) est à - 180 degrés.

    USE dimens_m, ONLY: iim
    use dynetat0_m, only: clon, grossismx, dzoomx, taux
    use invert_zoom_x_m, only: invert_zoom_x, nmax
    use nr_util, only: pi, pi_d, twopi, twopi_d, arth
    use principal_cshift_m, only: principal_cshift
    use tanh_cautious_m, only: tanh_cautious

    REAL, intent(out):: xprimm025(:), rlonv(:), xprimv(:) ! (iim + 1)
    real, intent(out):: rlonu(:), xprimu(:), xprimp025(:) ! (iim + 1)

    ! Local:
    real rlonm025(iim + 1), rlonp025(iim + 1), d_rlonv(iim)
    REAL dzoom, step
    DOUBLE PRECISION, dimension(0:nmax):: xtild, fhyp, G, Xf
    DOUBLE PRECISION ffdx, beta
    INTEGER i, is2
    DOUBLE PRECISION xxpr(nmax - 1), xmoy(nmax), fxm(nmax)

    !----------------------------------------------------------------------

    print *, "Call sequence information: fxhyp"

    test_grossismx: if (grossismx == 1.) then
       step = twopi / iim

       xprimm025(:iim) = step
       xprimp025(:iim) = step
       xprimv(:iim) = step
       xprimu(:iim) = step

       rlonv(:iim) = arth(- pi + clon, step, iim)
       rlonm025(:iim) = rlonv(:iim) - 0.25 * step
       rlonp025(:iim) = rlonv(:iim) + 0.25 * step
       rlonu(:iim) = rlonv(:iim) + 0.5 * step
    else test_grossismx
       dzoom = dzoomx * twopi_d
       xtild = arth(0d0, pi_d / nmax, nmax + 1)
       forall (i = 1:nmax) xmoy(i) = 0.5d0 * (xtild(i-1) + xtild(i))

       ! Compute fhyp:
       fhyp(1:nmax - 1) = tanh_cautious(taux * (dzoom / 2. &
            - xtild(1:nmax - 1)), xtild(1:nmax - 1) &
            * (pi_d - xtild(1:nmax - 1)))
       fhyp(0) = 1d0
       fhyp(nmax) = -1d0

       fxm = tanh_cautious(taux * (dzoom / 2. - xmoy), xmoy * (pi_d - xmoy))

       ! Calcul de beta 

       ffdx = 0.

       DO i = 1, nmax
          ffdx = ffdx + fxm(i) * (xtild(i) - xtild(i-1))
       END DO

       print *, "ffdx = ", ffdx
       beta = (pi_d - grossismx * ffdx) / (pi_d - ffdx)
       print *, "beta = ", beta

       IF (2. * beta - grossismx <= 0.) THEN
          print *, 'Bad choice of grossismx, taux, dzoomx.'
          print *, 'Decrease dzoomx or grossismx.'
          STOP 1
       END IF

       G = beta + (grossismx - beta) * fhyp

       ! Calcul de Xf

       xxpr = beta + (grossismx - beta) * fxm(:nmax - 1)
       Xf(0) = 0d0

       DO i = 1, nmax - 1
          Xf(i) = Xf(i-1) + xxpr(i) * (xtild(i) - xtild(i-1))
       END DO

       Xf(nmax) = pi_d

       call invert_zoom_x(xf, xtild, G, rlonm025(:iim), xprimm025(:iim), &
            xuv = - 0.25d0)
       call invert_zoom_x(xf, xtild, G, rlonv(:iim), xprimv(:iim), xuv = 0d0)
       call invert_zoom_x(xf, xtild, G, rlonu(:iim), xprimu(:iim), xuv = 0.5d0)
       call invert_zoom_x(xf, xtild, G, rlonp025(:iim), xprimp025(:iim), &
            xuv = 0.25d0)
    end if test_grossismx

    is2 = 0

    IF (MINval(rlonm025(:iim)) < - pi - 0.1 &
         .or. MAXval(rlonm025(:iim)) > pi + 0.1) THEN
       IF (clon <= 0.) THEN
          is2 = 1

          do while (rlonm025(is2) < - pi .and. is2 < iim)
             is2 = is2 + 1
          end do

          if (rlonm025(is2) < - pi) then
             print *, 'Rlonm025 plus petit que - pi !'
             STOP 1
          end if
       ELSE
          is2 = iim

          do while (rlonm025(is2) > pi .and. is2 > 1)
             is2 = is2 - 1
          end do

          if (rlonm025(is2) > pi) then
             print *, 'Rlonm025 plus grand que pi !'
             STOP 1
          end if
       END IF
    END IF

    call principal_cshift(is2, rlonm025, xprimm025)
    call principal_cshift(is2, rlonv, xprimv)
    call principal_cshift(is2, rlonu, xprimu)
    call principal_cshift(is2, rlonp025, xprimp025)

    forall (i = 1: iim) d_rlonv(i) = rlonv(i + 1) - rlonv(i)
    print *, "Minimum longitude step:", MINval(d_rlonv) * 180. / pi, "degrees"
    print *, "Maximum longitude step:", MAXval(d_rlonv) * 180. / pi, "degrees"

    ! Check that rlonm025 <= rlonv <= rlonp025 <= rlonu:
    DO i = 1, iim + 1
       IF (rlonp025(i) < rlonv(i)) THEN
          print *, 'rlonp025(', i, ') = ', rlonp025(i)
          print *, "< rlonv(", i, ") = ", rlonv(i)
          STOP 1
       END IF

       IF (rlonv(i) < rlonm025(i)) THEN 
          print *, 'rlonv(', i, ') = ', rlonv(i)
          print *, "< rlonm025(", i, ") = ", rlonm025(i)
          STOP 1
       END IF

       IF (rlonp025(i) > rlonu(i)) THEN
          print *, 'rlonp025(', i, ') = ', rlonp025(i)
          print *, "> rlonu(", i, ") = ", rlonu(i)
          STOP 1
       END IF
    END DO

  END SUBROUTINE fxhyp

end module fxhyp_m
