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
    ! 1., taux=0., clon=0.) est à - 180 degrés.

    USE dimens_m, ONLY: iim
    use fxhyp_loop_ik_m, only: fxhyp_loop_ik, nmax
    use nr_util, only: pi_d, twopi_d, arth
    use serre, only: clon, grossismx, dzoomx, taux

    REAL, intent(out):: xprimm025(:), rlonv(:), xprimv(:) ! (iim + 1)
    real, intent(out):: rlonu(:), xprimu(:), xprimp025(:) ! (iim + 1)

    ! Local:

    real rlonm025(iim + 1), rlonp025(iim + 1)
    REAL dzoom
    DOUBLE PRECISION xlon(iim)
    DOUBLE PRECISION xtild(0:2 * nmax)
    DOUBLE PRECISION fhyp(nmax:2 * nmax), ffdx, beta, Xprimt(0:2 * nmax)
    DOUBLE PRECISION Xf(0:2 * nmax), xxpr(2 * nmax)
    DOUBLE PRECISION xzoom, fa, fb
    INTEGER i
    DOUBLE PRECISION xmoy, fxm
    DOUBLE PRECISION decalx

    !----------------------------------------------------------------------

    print *, "Call sequence information: fxhyp"

    xzoom = clon * pi_d / 180d0

    IF (grossismx == 1.) THEN
       decalx = 1d0
    else
       decalx = 0.75d0
    END IF

    dzoom = dzoomx * twopi_d
    xtild = arth(- pi_d, pi_d / nmax, 2 * nmax + 1)

    ! Compute fhyp:
    DO i = nmax, 2 * nmax
       fa = taux * (dzoom / 2. - xtild(i))
       fb = xtild(i) * (pi_d - xtild(i))

       IF (200. * fb < - fa) THEN
          fhyp(i) = - 1.
       ELSE IF (200. * fb < fa) THEN
          fhyp(i) = 1.
       ELSE
          IF (ABS(fa) < 1e-13 .AND. ABS(fb) < 1e-13) THEN
             IF (200. * fb + fa < 1e-10) THEN
                fhyp(i) = - 1.
             ELSE IF (200. * fb - fa < 1e-10) THEN
                fhyp(i) = 1.
             END IF
          ELSE
             fhyp(i) = TANH(fa / fb)
          END IF
       END IF

       IF (xtild(i) == 0.) fhyp(i) = 1.
       IF (xtild(i) == pi_d) fhyp(i) = -1.
    END DO

    ! Calcul de beta 

    ffdx = 0.

    DO i = nmax + 1, 2 * nmax
       xmoy = 0.5 * (xtild(i-1) + xtild(i))
       fa = taux * (dzoom / 2. - xmoy)
       fb = xmoy * (pi_d - xmoy)

       IF (200. * fb < - fa) THEN
          fxm = - 1.
       ELSE IF (200. * fb < fa) THEN
          fxm = 1.
       ELSE
          IF (ABS(fa) < 1e-13 .AND. ABS(fb) < 1e-13) THEN
             IF (200. * fb + fa < 1e-10) THEN
                fxm = - 1.
             ELSE IF (200. * fb - fa < 1e-10) THEN
                fxm = 1.
             END IF
          ELSE
             fxm = TANH(fa / fb)
          END IF
       END IF

       IF (xmoy == 0.) fxm = 1.
       IF (xmoy == pi_d) fxm = -1.

       ffdx = ffdx + fxm * (xtild(i) - xtild(i-1))
    END DO

    print *, "ffdx = ", ffdx
    beta = (grossismx * ffdx - pi_d) / (ffdx - pi_d)
    print *, "beta = ", beta

    IF (2. * beta - grossismx <= 0.) THEN
       print *, 'Attention ! La valeur beta calculée dans fxhyp est mauvaise.'
       print *, 'Modifier les valeurs de grossismx, taux ou dzoomx et relancer.'
       STOP 1
    END IF

    ! calcul de Xprimt 
    Xprimt(nmax:2 * nmax) = beta + (grossismx - beta) * fhyp
    xprimt(:nmax - 1) = xprimt(2 * nmax:nmax + 1:- 1)

    ! Calcul de Xf

    DO i = nmax + 1, 2 * nmax
       xmoy = 0.5 * (xtild(i-1) + xtild(i))
       fa = taux * (dzoom / 2. - xmoy)
       fb = xmoy * (pi_d - xmoy)

       IF (200. * fb < - fa) THEN
          fxm = - 1.
       ELSE IF (200. * fb < fa) THEN
          fxm = 1.
       ELSE
          fxm = TANH(fa / fb)
       END IF

       IF (xmoy == 0.) fxm = 1.
       IF (xmoy == pi_d) fxm = -1.
       xxpr(i) = beta + (grossismx - beta) * fxm
    END DO

    xxpr(:nmax) = xxpr(2 * nmax:nmax + 1:- 1)

    Xf(0) = - pi_d

    DO i=1, 2 * nmax - 1
       Xf(i) = Xf(i-1) + xxpr(i) * (xtild(i) - xtild(i-1))
    END DO

    Xf(2 * nmax) = pi_d

    call fxhyp_loop_ik(1, decalx, xf, xtild, Xprimt, xzoom, rlonm025, &
         xprimm025, xuv = - 0.25d0)
    call fxhyp_loop_ik(2, decalx, xf, xtild, Xprimt, xzoom, rlonv, xprimv, &
         xuv = 0d0)
    call fxhyp_loop_ik(3, decalx, xf, xtild, Xprimt, xzoom, rlonu, xprimu, &
         xuv = 0.5d0)
    call fxhyp_loop_ik(4, decalx, xf, xtild, Xprimt, xzoom, rlonp025, &
         xprimp025, xuv = 0.25d0)

    print *

    forall (i = 1: iim) xlon(i) = rlonv(i + 1) - rlonv(i)
    print *, "Minimum longitude step:", MINval(xlon) * 180. / pi_d, "°"
    print *, "Maximum longitude step:", MAXval(xlon) * 180. / pi_d, "°"

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
