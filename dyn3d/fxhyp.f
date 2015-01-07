module fxhyp_m

  IMPLICIT NONE

contains

  SUBROUTINE fxhyp(xprimm025, rlonv, xprimv, rlonu, xprimu, xprimp025)

    ! From LMDZ4/libf/dyn3d/fxhyp.F, version 1.2, 2005/06/03 09:11:32
    ! Author: P. Le Van, from formulas by R. Sadourny

    ! Calcule les longitudes et dérivées dans la grille du GCM pour
    ! une fonction f(x) à dérivée tangente hyperbolique.

    ! On doit avoir grossismx \times dzoomx < pi (radians)

    USE dimens_m, ONLY: iim
    use nr_util, only: pi_d, twopi_d
    use serre, only: clon, grossismx, dzoomx, taux

    REAL, intent(out):: xprimm025(:), rlonv(:), xprimv(:) ! (iim + 1)
    real, intent(out):: rlonu(:), xprimu(:), xprimp025(:) ! (iim + 1)

    ! Local:

    DOUBLE PRECISION champmin, champmax
    real rlonm025(iim + 1), rlonp025(iim + 1)
    INTEGER, PARAMETER:: nmax = 30000, nmax2 = 2*nmax

    LOGICAL, PARAMETER:: scal180 = .TRUE.
    ! scal180 = .TRUE. si on veut avoir le premier point scalaire pour
    ! une grille reguliere (grossismx = 1., taux=0., clon=0.) a
    ! -180. degres. sinon scal180 = .FALSE.

    REAL dzoom
    DOUBLE PRECISION xlon(iim + 1), xprimm(iim + 1), xuv
    DOUBLE PRECISION xtild(0:nmax2)
    DOUBLE PRECISION fhyp(0:nmax2), ffdx, beta, Xprimt(0:nmax2)
    DOUBLE PRECISION Xf(0:nmax2), xxpr(0:nmax2)
    DOUBLE PRECISION xvrai(iim + 1), xxprim(iim + 1) 
    DOUBLE PRECISION my_eps, xzoom, fa, fb
    DOUBLE PRECISION Xf1, Xfi, a0, a1, a2, a3, xi2
    INTEGER i, it, ik, iter, ii, idif, ii1, ii2
    DOUBLE PRECISION xi, xo1, xmoy, xlon2, fxm, Xprimin
    DOUBLE PRECISION decalx
    INTEGER, save:: is2

    !----------------------------------------------------------------------

    my_eps = 1e-3
    xzoom = clon * pi_d / 180. 

    IF (grossismx == 1. .AND. scal180) THEN
       decalx = 1.
    else
       decalx = 0.75
    END IF

    IF (dzoomx < 1.) THEN
       dzoom = dzoomx * twopi_d
    ELSE IF (dzoomx < 25.) THEN
       print *, "Le paramètre dzoomx pour fxhyp est trop petit. " &
            // "L'augmenter et relancer."
       STOP 1
    ELSE
       dzoom = dzoomx * pi_d / 180.
    END IF

    print *, 'dzoom (rad):', dzoom

    DO i = 0, nmax2 
       xtild(i) = - pi_d + REAL(i) * twopi_d / nmax2
    END DO

    DO i = nmax, nmax2
       fa = taux* (dzoom / 2. - xtild(i))
       fb = xtild(i) * (pi_d - xtild(i))

       IF (200.* fb < - fa) THEN
          fhyp(i) = - 1.
       ELSE IF (200. * fb < fa) THEN
          fhyp(i) = 1.
       ELSE
          IF (ABS(fa) < 1e-13.AND.ABS(fb) < 1e-13) THEN
             IF (200.*fb + fa < 1e-10) THEN
                fhyp(i) = - 1.
             ELSE IF (200.*fb - fa < 1e-10) THEN
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

    DO i = nmax + 1, nmax2
       xmoy = 0.5 * (xtild(i-1) + xtild(i))
       fa = taux* (dzoom / 2. - xmoy)
       fb = xmoy * (pi_d - xmoy)

       IF (200.* fb < - fa) THEN
          fxm = - 1.
       ELSE IF (200. * fb < fa) THEN
          fxm = 1.
       ELSE
          IF (ABS(fa) < 1e-13.AND.ABS(fb) < 1e-13) THEN
             IF (200.*fb + fa < 1e-10) THEN
                fxm = - 1.
             ELSE IF (200.*fb - fa < 1e-10) THEN
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

    beta = (grossismx * ffdx - pi_d) / (ffdx - pi_d)

    IF (2. * beta - grossismx <= 0.) THEN
       print *, 'Attention ! La valeur beta calculée dans fxhyp est mauvaise.'
       print *, 'Modifier les valeurs de grossismx, taux ou dzoomx et relancer.'
       STOP 1
    END IF

    ! calcul de Xprimt 

    DO i = nmax, nmax2
       Xprimt(i) = beta + (grossismx - beta) * fhyp(i)
    END DO

    DO i = nmax + 1, nmax2
       Xprimt(nmax2 - i) = Xprimt(i)
    END DO

    ! Calcul de Xf 

    Xf(0) = - pi_d

    DO i = nmax + 1, nmax2
       xmoy = 0.5 * (xtild(i-1) + xtild(i))
       fa = taux* (dzoom / 2. - xmoy)
       fb = xmoy * (pi_d - xmoy)

       IF (200.* fb < - fa) THEN
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

    DO i = nmax + 1, nmax2
       xxpr(nmax2-i + 1) = xxpr(i)
    END DO

    DO i=1, nmax2
       Xf(i) = Xf(i-1) + xxpr(i) * (xtild(i) - xtild(i-1))
    END DO

    ! xuv = 0. si calcul aux points scalaires 
    ! xuv = 0.5 si calcul aux points U 

    loop_ik: DO ik = 1, 4
       IF (ik == 1) THEN
          xuv = -0.25
       ELSE IF (ik == 2) THEN
          xuv = 0.
       ELSE IF (ik == 3) THEN
          xuv = 0.50
       ELSE IF (ik == 4) THEN
          xuv = 0.25
       END IF

       xo1 = 0.

       ii1=1
       ii2=iim
       IF (ik == 1.and.grossismx == 1.) THEN
          ii1 = 2 
          ii2 = iim + 1
       END IF

       DO i = ii1, ii2
          xlon2 = - pi_d + (REAL(i) + xuv - decalx) * twopi_d / REAL(iim) 
          Xfi = xlon2

          it = nmax2
          do while (xfi < xf(it) .and. it >= 1)
             it = it - 1
          end do

          ! Calcul de Xf(xi) 

          xi = xtild(it)

          IF (it == nmax2) THEN
             it = nmax2 -1
             Xf(it + 1) = pi_d
          END IF

          ! Appel de la routine qui calcule les coefficients a0, a1,
          ! a2, a3 d'un polynome de degre 3 qui passe par les points
          ! (Xf(it), xtild(it)) et (Xf(it + 1), xtild(it + 1))

          CALL coefpoly(Xf(it), Xf(it + 1), Xprimt(it), Xprimt(it + 1), &
               xtild(it), xtild(it + 1), a0, a1, a2, a3)

          Xf1 = Xf(it)
          Xprimin = a1 + 2.* a2 * xi + 3.*a3 * xi *xi

          iter = 1

          do
             xi = xi - (Xf1 - Xfi) / Xprimin
             IF (ABS(xi - xo1) <= my_eps .or. iter == 300) exit
             xo1 = xi
             xi2 = xi * xi
             Xf1 = a0 + a1 * xi + a2 * xi2 + a3 * xi2 * xi
             Xprimin = a1 + 2.* a2 * xi + 3.* a3 * xi2
          end DO

          if (ABS(xi - xo1) > my_eps) then
             ! iter == 300
             print *, 'Pas de solution.'
             print *, i, xlon2
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
             print *, 'Problème avec rlonu(', i + 1, &
                  ') plus petit que rlonu(', i, ')'
             STOP 1
          END IF
       END DO

       ! Réorganisation des longitudes pour les avoir entre - pi et pi 

       champmin = 1e12
       champmax = -1e12
       DO i = 1, iim
          champmin = MIN(champmin, xvrai(i))
          champmax = MAX(champmax, xvrai(i))
       END DO

       IF (.not. (champmin >= -pi_d - 0.1 .and. champmax <= pi_d + 0.1)) THEN
          print *, 'Reorganisation des longitudes pour avoir entre - pi', &
               ' et pi '

          IF (xzoom <= 0.) THEN
             IF (ik == 1) THEN
                i = 1

                do while (xvrai(i) < - pi_d .and. i < iim)
                   i = i + 1
                end do

                if (xvrai(i) < - pi_d) then
                   print *, 'Xvrai plus petit que - pi !'
                   STOP 1
                end if

                is2 = i
             END IF

             IF (is2.NE. 1) THEN
                DO ii = is2, iim
                   xlon(ii-is2 + 1) = xvrai(ii)
                   xprimm(ii-is2 + 1) = xxprim(ii)
                END DO
                DO ii = 1, is2 -1
                   xlon(ii + iim-is2 + 1) = xvrai(ii) + twopi_d
                   xprimm(ii + iim-is2 + 1) = xxprim(ii) 
                END DO
             END IF
          ELSE 
             IF (ik == 1) THEN
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

             idif = iim -is2

             DO ii = 1, is2
                xlon(ii + idif) = xvrai(ii)
                xprimm(ii + idif) = xxprim(ii)
             END DO

             DO ii = 1, idif
                xlon(ii) = xvrai(ii + is2) - twopi_d
                xprimm(ii) = xxprim(ii + is2) 
             END DO
          END IF
       END IF

       ! Fin de la reorganisation 

       xlon(iim + 1) = xlon(1) + twopi_d
       xprimm(iim + 1) = xprimm(1)

       DO i = 1, iim + 1
          xvrai(i) = xlon(i)*180. / pi_d
       END DO

       IF (ik == 1) THEN
          DO i = 1, iim + 1
             rlonm025(i) = xlon(i)
             xprimm025(i) = xprimm(i)
          END DO
       ELSE IF (ik == 2) THEN
          rlonv = xlon
          xprimv = xprimm
       ELSE IF (ik == 3) THEN
          DO i = 1, iim + 1
             rlonu(i) = xlon(i)
             xprimu(i) = xprimm(i)
          END DO
       ELSE IF (ik == 4) THEN
          DO i = 1, iim + 1
             rlonp025(i) = xlon(i)
             xprimp025(i) = xprimm(i)
          END DO
       END IF
    end DO loop_ik

    print *

    DO i = 1, iim
       xlon(i) = rlonv(i + 1) - rlonv(i)
    END DO
    champmin = 1e12
    champmax = -1e12
    DO i = 1, iim
       champmin = MIN(champmin, xlon(i))
       champmax = MAX(champmax, xlon(i))
    END DO
    champmin = champmin * 180. / pi_d
    champmax = champmax * 180. / pi_d

    DO i = 1, iim + 1
       IF (rlonp025(i) < rlonv(i)) THEN
          print *, ' Attention ! rlonp025 < rlonv', i
          STOP 1
       END IF

       IF (rlonv(i) < rlonm025(i)) THEN 
          print *, ' Attention ! rlonm025 > rlonv', i
          STOP 1
       END IF

       IF (rlonp025(i) > rlonu(i)) THEN
          print *, ' Attention ! rlonp025 > rlonu', i
          STOP 1
       END IF
    END DO

    print *, ' Longitudes '
    print 3, champmin, champmax

3   Format(1x, ' Au centre du zoom, la longueur de la maille est', &
         ' d environ ', f0.2, ' degres ', /, &
         ' alors que la maille en dehors de la zone du zoom est ', &
         "d'environ", f0.2, ' degres ')

  END SUBROUTINE fxhyp

end module fxhyp_m
