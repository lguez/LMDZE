module fyhyp_m

  IMPLICIT NONE

contains

  SUBROUTINE fyhyp(rlatu, yyprimu, rlatv, rlatu2, yprimu2, rlatu1, yprimu1)

    ! From LMDZ4/libf/dyn3d/fyhyp.F, version 1.2, 2005/06/03 09:11:32

    ! Author: P. Le Van, from analysis by R. Sadourny 

    ! Calcule les latitudes et dérivées dans la grille du GCM pour une
    ! fonction f(y) à dérivée tangente hyperbolique.

    ! Il vaut mieux avoir : grossismy * dzoom < pi / 2

    use coefpoly_m, only: coefpoly
    USE dimens_m, only: jjm
    use heavyside_m, only: heavyside
    use serre, only: clat, grossismy, dzoomy, tauy

    REAL, intent(out):: rlatu(jjm + 1), yyprimu(jjm + 1)
    REAL, intent(out):: rlatv(jjm)
    real, intent(out):: rlatu2(jjm), yprimu2(jjm), rlatu1(jjm), yprimu1(jjm)

    ! Local: 

    DOUBLE PRECISION champmin, champmax
    INTEGER, PARAMETER:: nmax=30000, nmax2=2*nmax
    REAL dzoom ! distance totale de la zone du zoom (en radians)
    DOUBLE PRECISION ylat(jjm + 1), yprim(jjm + 1)
    DOUBLE PRECISION yuv
    DOUBLE PRECISION, save:: yt(0:nmax2)
    DOUBLE PRECISION fhyp(0:nmax2), beta
    DOUBLE PRECISION, save:: ytprim(0:nmax2)
    DOUBLE PRECISION fxm(0:nmax2)
    DOUBLE PRECISION, save:: yf(0:nmax2)
    DOUBLE PRECISION yypr(0:nmax2)
    DOUBLE PRECISION yvrai(jjm + 1), yprimm(jjm + 1), ylatt(jjm + 1)
    DOUBLE PRECISION pi, pis2, epsilon, pisjm
    DOUBLE PRECISION yo1, yi, ylon2, ymoy, yprimin
    DOUBLE PRECISION yfi, yf1, ffdy
    DOUBLE PRECISION ypn, deply, y00
    SAVE y00, deply

    INTEGER i, j, it, ik, iter, jlat
    INTEGER jpn, jjpn
    SAVE jpn
    DOUBLE PRECISION a0, a1, a2, a3, yi2, heavyy0, heavyy0m
    DOUBLE PRECISION fa(0:nmax2), fb(0:nmax2)
    REAL y0min, y0max

    !-------------------------------------------------------------------

    print *, "Call sequence information: fyhyp"

    pi = 2.*asin(1.)
    pis2 = pi/2.
    pisjm = pi/real(jjm)
    epsilon = 1e-3
    dzoom = dzoomy*pi
    print *, 'yzoom(rad), grossismy, tauy, dzoom (rad):'
    print *, clat, grossismy, tauy, dzoom

    DO i = 0, nmax2
       yt(i) = -pis2 + real(i)*pi/nmax2
    END DO

    heavyy0m = heavyside(-clat)
    heavyy0 = heavyside(clat)
    y0min = 2.*clat*heavyy0m - pis2
    y0max = 2.*clat*heavyy0 + pis2

    fa = 999.999
    fb = 999.999

    DO i = 0, nmax2
       IF (yt(i)<clat) THEN
          fa(i) = tauy*(yt(i)-clat + dzoom/2.)
          fb(i) = (yt(i)-2.*clat*heavyy0m + pis2)*(clat-yt(i))
       ELSE IF (yt(i)>clat) THEN
          fa(i) = tauy*(clat-yt(i) + dzoom/2.)
          fb(i) = (2.*clat*heavyy0-yt(i) + pis2)*(yt(i)-clat)
       END IF

       IF (200.*fb(i)<-fa(i)) THEN
          fhyp(i) = -1.
       ELSE IF (200.*fb(i)<fa(i)) THEN
          fhyp(i) = 1.
       ELSE
          fhyp(i) = tanh(fa(i)/fb(i))
       END IF

       IF (yt(i)==clat) fhyp(i) = 1.
       IF (yt(i)==y0min .OR. yt(i)==y0max) fhyp(i) = -1.
    END DO

    ! Calcul de beta 

    ffdy = 0.

    DO i = 1, nmax2
       ymoy = 0.5*(yt(i-1) + yt(i))
       IF (ymoy<clat) THEN
          fa(i) = tauy*(ymoy-clat + dzoom/2.)
          fb(i) = (ymoy-2.*clat*heavyy0m + pis2)*(clat-ymoy)
       ELSE IF (ymoy>clat) THEN
          fa(i) = tauy*(clat-ymoy + dzoom/2.)
          fb(i) = (2.*clat*heavyy0-ymoy + pis2)*(ymoy-clat)
       END IF

       IF (200.*fb(i)<-fa(i)) THEN
          fxm(i) = -1.
       ELSE IF (200.*fb(i)<fa(i)) THEN
          fxm(i) = 1.
       ELSE
          fxm(i) = tanh(fa(i)/fb(i))
       END IF
       IF (ymoy==clat) fxm(i) = 1.
       IF (ymoy==y0min .OR. yt(i)==y0max) fxm(i) = -1.
       ffdy = ffdy + fxm(i)*(yt(i)-yt(i-1))
    END DO

    beta = (grossismy*ffdy-pi)/(ffdy-pi)

    IF (2. * beta - grossismy <= 0.) THEN
       print *, 'Attention ! La valeur beta calculee dans la routine fyhyp ' &
            // 'est mauvaise. Modifier les valeurs de grossismy, tauy ou ' &
            // 'dzoomy et relancer.'
       STOP 1
    END IF

    ! calcul de Ytprim 

    DO i = 0, nmax2
       ytprim(i) = beta + (grossismy-beta)*fhyp(i)
    END DO

    ! Calcul de Yf 

    yf(0) = -pis2
    DO i = 1, nmax2
       yypr(i) = beta + (grossismy-beta)*fxm(i)
    END DO

    DO i = 1, nmax2
       yf(i) = yf(i-1) + yypr(i)*(yt(i)-yt(i-1))
    END DO

    ! yuv = 0. si calcul des latitudes aux pts. U 
    ! yuv = 0.5 si calcul des latitudes aux pts. V 

    loop_ik: DO ik = 1, 4
       IF (ik==1) THEN
          yuv = 0.
          jlat = jjm + 1
       ELSE IF (ik==2) THEN
          yuv = 0.5
          jlat = jjm
       ELSE IF (ik==3) THEN
          yuv = 0.25
          jlat = jjm
       ELSE IF (ik==4) THEN
          yuv = 0.75
          jlat = jjm
       END IF

       yo1 = 0.
       DO j = 1, jlat
          yo1 = 0.
          ylon2 = -pis2 + pisjm*(real(j) + yuv-1.)
          yfi = ylon2

          it = nmax2
          DO while (it >= 1 .and. yfi < yf(it))
             it = it - 1
          END DO

          yi = yt(it)
          IF (it==nmax2) THEN
             it = nmax2 - 1
             yf(it + 1) = pis2
          END IF

          ! Interpolation entre yi(it) et yi(it + 1) pour avoir Y(yi)
          ! et Y'(yi) 

          CALL coefpoly(yf(it), yf(it + 1), ytprim(it), ytprim(it + 1), &
               yt(it), yt(it + 1), a0, a1, a2, a3)

          yf1 = yf(it)
          yprimin = a1 + 2.*a2*yi + 3.*a3*yi*yi

          iter = 1
          DO
             yi = yi - (yf1-yfi)/yprimin
             IF (abs(yi-yo1)<=epsilon .or. iter == 300) exit
             yo1 = yi
             yi2 = yi*yi
             yf1 = a0 + a1*yi + a2*yi2 + a3*yi2*yi
             yprimin = a1 + 2.*a2*yi + 3.*a3*yi2
          END DO
          if (abs(yi-yo1) > epsilon) then
             print *, 'Pas de solution.', j, ylon2
             STOP 1
          end if

          yprimin = a1 + 2.*a2*yi + 3.*a3*yi*yi
          yprim(j) = pi/(jjm*yprimin)
          yvrai(j) = yi
       END DO

       DO j = 1, jlat - 1
          IF (yvrai(j + 1)<yvrai(j)) THEN
             print *, 'Problème avec rlat(', j + 1, ') plus petit que rlat(', &
                  j, ')'
             STOP 1
          END IF
       END DO

       print *, 'Reorganisation des latitudes pour avoir entre - pi/2 et pi/2'

       IF (ik==1) THEN
          ypn = pis2
          DO j = jjm + 1, 1, -1
             IF (yvrai(j)<=ypn) exit
          END DO

          jpn = j
          y00 = yvrai(jpn)
          deply = pis2 - y00
       END IF

       DO j = 1, jjm + 1 - jpn
          ylatt(j) = -pis2 - y00 + yvrai(jpn + j-1)
          yprimm(j) = yprim(jpn + j-1)
       END DO

       jjpn = jpn
       IF (jlat==jjm) jjpn = jpn - 1

       DO j = 1, jjpn
          ylatt(j + jjm + 1-jpn) = yvrai(j) + deply
          yprimm(j + jjm + 1-jpn) = yprim(j)
       END DO

       ! Fin de la reorganisation

       DO j = 1, jlat
          ylat(j) = ylatt(jlat + 1-j)
          yprim(j) = yprimm(jlat + 1-j)
       END DO

       DO j = 1, jlat
          yvrai(j) = ylat(j)*180./pi
       END DO

       IF (ik==1) THEN
          DO j = 1, jjm + 1
             rlatu(j) = ylat(j)
             yyprimu(j) = yprim(j)
          END DO
       ELSE IF (ik==2) THEN
          DO j = 1, jjm
             rlatv(j) = ylat(j)
          END DO
       ELSE IF (ik==3) THEN
          DO j = 1, jjm
             rlatu2(j) = ylat(j)
             yprimu2(j) = yprim(j)
          END DO
       ELSE IF (ik==4) THEN
          DO j = 1, jjm
             rlatu1(j) = ylat(j)
             yprimu1(j) = yprim(j)
          END DO
       END IF
    END DO loop_ik

    DO j = 1, jjm
       ylat(j) = rlatu(j) - rlatu(j + 1)
    END DO
    champmin = 1e12
    champmax = -1e12
    DO j = 1, jjm
       champmin = min(champmin, ylat(j))
       champmax = max(champmax, ylat(j))
    END DO
    champmin = champmin*180./pi
    champmax = champmax*180./pi

    DO j = 1, jjm
       IF (rlatu1(j) <= rlatu2(j)) THEN
          print *, 'Attention ! rlatu1 < rlatu2 ', rlatu1(j), rlatu2(j), j
          STOP 13
       ENDIF

       IF (rlatu2(j) <= rlatu(j+1)) THEN
          print *, 'Attention ! rlatu2 < rlatup1 ', rlatu2(j), rlatu(j+1), j
          STOP 14
       ENDIF

       IF (rlatu(j) <= rlatu1(j)) THEN
          print *, ' Attention ! rlatu < rlatu1 ', rlatu(j), rlatu1(j), j
          STOP 15
       ENDIF

       IF (rlatv(j) <= rlatu2(j)) THEN
          print *, ' Attention ! rlatv < rlatu2 ', rlatv(j), rlatu2(j), j
          STOP 16
       ENDIF

       IF (rlatv(j) >= rlatu1(j)) THEN
          print *, ' Attention ! rlatv > rlatu1 ', rlatv(j), rlatu1(j), j
          STOP 17
       ENDIF

       IF (rlatv(j) >= rlatu(j)) THEN
          print *, ' Attention ! rlatv > rlatu ', rlatv(j), rlatu(j), j
          STOP 18
       ENDIF
    ENDDO

    print *, 'Latitudes'
    print 3, champmin, champmax

3   Format(1x, ' Au centre du zoom, la longueur de la maille est', &
         ' d environ ', f0.2, ' degres ', /, &
         ' alors que la maille en dehors de la zone du zoom est ', &
         "d'environ ", f0.2, ' degres ')

  END SUBROUTINE fyhyp

end module fyhyp_m
