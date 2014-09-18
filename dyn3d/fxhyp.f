module fxhyp_m

  IMPLICIT NONE

contains

  SUBROUTINE fxhyp(xzoomdeg, grossism, dzooma, tau, rlonm025, xprimm025, &
       rlonv, xprimv, rlonu, xprimu, rlonp025, xprimp025, champmin, champmax)

    ! From LMDZ4/libf/dyn3d/fxhyp.F, version 1.2, 2005/06/03 09:11:32
    ! Author: P. Le Van 

    ! Calcule les longitudes et dérivées dans la grille du GCM pour
    ! une fonction f(x) à tangente hyperbolique.

    ! On doit avoir grossism \times dzoom < pi (radians), en longitude.

    USE dimens_m, ONLY: iim
    USE paramet_m, ONLY: iip1

    REAL, intent(in):: xzoomdeg

    REAL, intent(in):: grossism
    ! grossissement (= 2 si 2 fois, = 3 si 3 fois, etc.)

    REAL, intent(in):: dzooma ! distance totale de la zone du zoom

    REAL, intent(in):: tau
    ! raideur de la transition de l'intérieur à l'extérieur du zoom

    ! arguments de sortie 

    REAL, dimension(iip1):: rlonm025, xprimm025, rlonv, xprimv
    real, dimension(iip1):: rlonu, xprimu, rlonp025, xprimp025

    DOUBLE PRECISION, intent(out):: champmin, champmax

    ! Local:

    INTEGER, PARAMETER:: nmax = 30000, nmax2 = 2*nmax

    LOGICAL, PARAMETER:: scal180 = .TRUE.
    ! scal180 = .TRUE. si on veut avoir le premier point scalaire pour
    ! une grille reguliere (grossism = 1., tau=0., clon=0.) a
    ! -180. degres. sinon scal180 = .FALSE.

    REAL dzoom
    DOUBLE PRECISION xlon(iip1), xprimm(iip1), xuv
    DOUBLE PRECISION xtild(0:nmax2)
    DOUBLE PRECISION fhyp(0:nmax2), ffdx, beta, Xprimt(0:nmax2)
    DOUBLE PRECISION Xf(0:nmax2), xxpr(0:nmax2)
    DOUBLE PRECISION xvrai(iip1), xxprim(iip1) 
    DOUBLE PRECISION pi, depi, epsilon, xzoom, fa, fb
    DOUBLE PRECISION Xf1, Xfi, a0, a1, a2, a3, xi2
    INTEGER i, it, ik, iter, ii, idif, ii1, ii2
    DOUBLE PRECISION xi, xo1, xmoy, xlon2, fxm, Xprimin
    DOUBLE PRECISION decalx
    INTEGER is2
    SAVE is2

    !----------------------------------------------------------------------

    pi = 2. * ASIN(1.)
    depi = 2. * pi
    epsilon = 1.e-3
    xzoom = xzoomdeg * pi/180. 

    decalx = .75
    IF (grossism == 1. .AND. scal180) THEN
       decalx = 1.
    ENDIF

    print *, 'FXHYP scal180, decalx', scal180, decalx

    IF (dzooma.LT.1.) THEN
       dzoom = dzooma * depi
    ELSEIF (dzooma.LT. 25.) THEN
       print *, "Le paramètre dzoomx pour fxhyp est trop petit. " &
            // "L'augmenter et relancer."
       STOP 1
    ELSE
       dzoom = dzooma * pi/180.
    END IF

    print *, ' xzoom(rad), grossism, tau, dzoom (rad):'
    print *, xzoom, grossism, tau, dzoom

    DO i = 0, nmax2 
       xtild(i) = - pi + REAL(i) * depi /nmax2
    ENDDO

    DO i = nmax, nmax2
       fa = tau* (dzoom/2. - xtild(i))
       fb = xtild(i) * (pi - xtild(i))

       IF (200.* fb .LT. - fa) THEN
          fhyp (i) = - 1.
       ELSEIF (200. * fb .LT. fa) THEN
          fhyp (i) = 1.
       ELSE
          IF (ABS(fa).LT.1.e-13.AND.ABS(fb).LT.1.e-13) THEN
             IF (200.*fb + fa.LT.1.e-10) THEN
                fhyp (i) = - 1.
             ELSEIF (200.*fb - fa.LT.1.e-10) THEN
                fhyp (i) = 1.
             ENDIF
          ELSE
             fhyp (i) = TANH (fa/fb)
          ENDIF
       END IF

       IF (xtild(i) == 0.) fhyp(i) = 1.
       IF (xtild(i) == pi) fhyp(i) = -1.
    END DO

    ! Calcul de beta 

    ffdx = 0.

    DO i = nmax + 1, nmax2
       xmoy = 0.5 * (xtild(i-1) + xtild(i))
       fa = tau* (dzoom/2. - xmoy)
       fb = xmoy * (pi - xmoy)

       IF (200.* fb .LT. - fa) THEN
          fxm = - 1.
       ELSEIF (200. * fb .LT. fa) THEN
          fxm = 1.
       ELSE
          IF (ABS(fa).LT.1.e-13.AND.ABS(fb).LT.1.e-13) THEN
             IF (200.*fb + fa.LT.1.e-10) THEN
                fxm = - 1.
             ELSEIF (200.*fb - fa.LT.1.e-10) THEN
                fxm = 1.
             ENDIF
          ELSE
             fxm = TANH (fa/fb)
          ENDIF
       ENDIF

       IF (xmoy == 0.) fxm = 1.
       IF (xmoy == pi) fxm = -1.

       ffdx = ffdx + fxm * (xtild(i) - xtild(i-1))
    ENDDO

    beta = (grossism * ffdx - pi) / (ffdx - pi)

    IF (2.*beta - grossism <= 0.) THEN
       print *, 'Attention ! La valeur beta calculée dans fxhyp est mauvaise.'
       print *, 'Modifier les valeurs de grossismx, tau ou dzoomx et relancer.'
       STOP 1
    END IF

    ! calcul de Xprimt 

    DO i = nmax, nmax2
       Xprimt(i) = beta + (grossism - beta) * fhyp(i)
    END DO

    DO i = nmax + 1, nmax2
       Xprimt(nmax2 - i) = Xprimt(i)
    END DO

    ! Calcul de Xf 

    Xf(0) = - pi

    DO i = nmax + 1, nmax2
       xmoy = 0.5 * (xtild(i-1) + xtild(i))
       fa = tau* (dzoom/2. - xmoy)
       fb = xmoy * (pi - xmoy)

       IF (200.* fb .LT. - fa) THEN
          fxm = - 1.
       ELSEIF (200. * fb .LT. fa) THEN
          fxm = 1.
       ELSE
          fxm = TANH (fa/fb)
       ENDIF

       IF (xmoy == 0.) fxm = 1.
       IF (xmoy == pi) fxm = -1.
       xxpr(i) = beta + (grossism - beta) * fxm
    ENDDO

    DO i = nmax + 1, nmax2
       xxpr(nmax2-i + 1) = xxpr(i)
    ENDDO

    DO i=1, nmax2
       Xf(i) = Xf(i-1) + xxpr(i) * (xtild(i) - xtild(i-1))
    ENDDO

    ! xuv = 0. si calcul aux pts scalaires 
    ! xuv = 0.5 si calcul aux pts U 

    print *

    DO ik = 1, 4
       IF (ik == 1) THEN
          xuv = -0.25
       ELSE IF (ik == 2) THEN
          xuv = 0.
       ELSE IF (ik == 3) THEN
          xuv = 0.50
       ELSE IF (ik == 4) THEN
          xuv = 0.25
       ENDIF

       xo1 = 0.

       ii1=1
       ii2=iim
       IF (ik == 1.and.grossism == 1.) THEN
          ii1 = 2 
          ii2 = iim + 1
       ENDIF

       DO i = ii1, ii2
          xlon2 = - pi + (REAL(i) + xuv - decalx) * depi / REAL(iim) 
          Xfi = xlon2

          it = nmax2
          do while (xfi < xf(it) .and. it >= 1)
             it = it - 1
          end do

          ! Calcul de Xf(xi) 

          xi = xtild(it)

          IF (it == nmax2) THEN
             it = nmax2 -1
             Xf(it + 1) = pi
          ENDIF

          ! Appel de la routine qui calcule les coefficients a0, a1,
          ! a2, a3 d'un polynome de degre 3 qui passe par les points
          ! (Xf(it), xtild(it)) et (Xf(it + 1), xtild(it + 1))

          CALL coefpoly(Xf(it), Xf(it + 1), Xprimt(it), Xprimt(it + 1), &
               xtild(it), xtild(it + 1), a0, a1, a2, a3)

          Xf1 = Xf(it)
          Xprimin = a1 + 2.* a2 * xi + 3.*a3 * xi *xi

          iter = 1

          do
             xi = xi - (Xf1 - Xfi)/ Xprimin
             IF (ABS(xi - xo1) <= epsilon .or. iter == 300) exit
             xo1 = xi
             xi2 = xi * xi
             Xf1 = a0 + a1 * xi + a2 * xi2 + a3 * xi2 * xi
             Xprimin = a1 + 2.* a2 * xi + 3.* a3 * xi2
          end DO

          if (ABS(xi - xo1) > epsilon) then
             ! iter == 300
             print *, 'Pas de solution.'
             print *, i, xlon2
             STOP 1
          end if


          xxprim(i) = depi/ (REAL(iim) * Xprimin)
          xvrai(i) = xi + xzoom
       end DO

       IF (ik == 1.and.grossism == 1.) THEN
          xvrai(1) = xvrai(iip1)-depi
          xxprim(1) = xxprim(iip1)
       ENDIF
       DO i = 1, iim
          xlon(i) = xvrai(i)
          xprimm(i) = xxprim(i)
       ENDDO
       DO i = 1, iim -1
          IF (xvrai(i + 1).LT. xvrai(i)) THEN
             print *, 'Problème avec rlonu(', i + 1, &
                  ') plus petit que rlonu(', i, ')'
             STOP 1
          ENDIF
       ENDDO

       ! Reorganisation des longitudes pour les avoir entre - pi et pi 

       champmin = 1.e12
       champmax = -1.e12
       DO i = 1, iim
          champmin = MIN(champmin, xvrai(i))
          champmax = MAX(champmax, xvrai(i))
       ENDDO

       IF (.not. (champmin >= -pi-0.10.and.champmax <= pi + 0.10)) THEN
          print *, 'Reorganisation des longitudes pour avoir entre - pi', &
               ' et pi '

          IF (xzoom <= 0.) THEN
             IF (ik == 1) THEN
                i = 1

                do while (xvrai(i) < - pi .and. i < iim)
                   i = i + 1
                end do

                if (xvrai(i) < - pi) then
                   print *, ' PBS. 1 ! Xvrai plus petit que - pi ! '
                   STOP 1
                end if

                is2 = i
             ENDIF

             IF (is2.NE. 1) THEN
                DO ii = is2, iim
                   xlon (ii-is2 + 1) = xvrai(ii)
                   xprimm(ii-is2 + 1) = xxprim(ii)
                ENDDO
                DO ii = 1, is2 -1
                   xlon (ii + iim-is2 + 1) = xvrai(ii) + depi
                   xprimm(ii + iim-is2 + 1) = xxprim(ii) 
                ENDDO
             ENDIF
          ELSE 
             IF (ik == 1) THEN
                i = iim

                do while (xvrai(i) > pi .and. i > 1)
                   i = i - 1
                end do

                if (xvrai(i) > pi) then
                   print *, ' PBS. 2 ! Xvrai plus grand que pi ! '
                   STOP 1
                end if

                is2 = i
             ENDIF
             idif = iim -is2
             DO ii = 1, is2
                xlon (ii + idif) = xvrai(ii)
                xprimm(ii + idif) = xxprim(ii)
             ENDDO
             DO ii = 1, idif
                xlon (ii) = xvrai (ii + is2) - depi
                xprimm(ii) = xxprim(ii + is2) 
             ENDDO
          ENDIF
       ENDIF

       ! Fin de la reorganisation 

       xlon (iip1) = xlon(1) + depi
       xprimm(iip1) = xprimm (1)

       DO i = 1, iim + 1
          xvrai(i) = xlon(i)*180./pi
       ENDDO

       IF (ik == 1) THEN
          DO i = 1, iim + 1
             rlonm025(i) = xlon(i)
             xprimm025(i) = xprimm(i)
          ENDDO
       ELSE IF (ik == 2) THEN
          DO i = 1, iim + 1
             rlonv(i) = xlon(i)
             xprimv(i) = xprimm(i)
          ENDDO
       ELSE IF (ik == 3) THEN
          DO i = 1, iim + 1
             rlonu(i) = xlon(i)
             xprimu(i) = xprimm(i)
          ENDDO
       ELSE IF (ik == 4) THEN
          DO i = 1, iim + 1
             rlonp025(i) = xlon(i)
             xprimp025(i) = xprimm(i)
          ENDDO
       ENDIF
    end DO

    print *

    DO i = 1, iim
       xlon(i) = rlonv(i + 1) - rlonv(i)
    ENDDO
    champmin = 1.e12
    champmax = -1.e12
    DO i = 1, iim
       champmin = MIN(champmin, xlon(i))
       champmax = MAX(champmax, xlon(i))
    ENDDO
    champmin = champmin * 180./pi
    champmax = champmax * 180./pi

  END SUBROUTINE fxhyp

end module fxhyp_m
