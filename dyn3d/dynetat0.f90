module dynetat0_m

  IMPLICIT NONE

  INTEGER, protected, save:: day_ini 
  ! day number at the beginning of the run, based at value 1 on
  ! January 1st of annee_ref

  real, protected, save, allocatable:: rlatu(:) ! (jjm + 1)
  ! latitudes of points of the "scalar" and "u" grid, in rad

  real, protected, save, allocatable:: rlatv(:) ! (jjm) 
  ! latitudes of points of the "v" grid, in rad, in decreasing order

  real, protected, save, allocatable:: rlonu(:) ! (iim + 1)
  ! longitudes of points of the "u" grid, in rad

  real, protected, save, allocatable:: rlonv(:) ! (iim + 1)
  ! longitudes of points of the "scalar" and "v" grid, in rad

  real, protected, save, allocatable:: xprimu(:), xprimv(:) ! (iim + 1)
  ! 2 pi / iim * (derivative of the longitudinal zoom function)(rlon[uv])

  REAL, protected, save, allocatable:: xprimm025(:), xprimp025(:) ! (iim + 1)
  REAL, protected, save, allocatable:: rlatu1(:), rlatu2(:), yprimu1(:), yprimu2(:) ! (jjm)
  REAL, save:: ang0, etot0, ptot0, ztot0, stot0
  INTEGER, PARAMETER, private:: nmax = 30000
  INTEGER, protected, save:: itau_dyn

contains

  SUBROUTINE dynetat0(vcov, ucov, teta, q, masse, ps)

    ! From dynetat0.F, version 1.2, 2004/06/22 11:45:30
    ! Authors: P. Le Van, L. Fairhead
    ! This procedure reads the initial state of the atmosphere.

    ! Libraries:
    use netcdf, only: NF90_NOWRITE, NF90_NOERR
    use netcdf95, only: NF95_GET_VAR, nf95_open, nf95_inq_varid, NF95_CLOSE, &
         NF95_Gw_VAR
    use jumble, only: assert

    use conf_gcm_m, only: raz_date
    use dimensions, only: iim, jjm, llm, nqmx
    use dynetat0_chosen_m, only: day_ref
    use infotrac_init_m, only: tname

    REAL, intent(out):: vcov(: , :, :) ! (iim + 1, jjm, llm)
    REAL, intent(out):: ucov(:, :, :) ! (iim + 1, jjm + 1, llm)
    REAL, intent(out):: teta(:, :, :) ! (iim + 1, jjm + 1, llm)
    REAL, intent(out):: q(:, :, :, :) ! (iim + 1, jjm + 1, llm, nqmx)
    REAL, intent(out):: masse(:, :, :) ! (iim + 1, jjm + 1, llm)
    REAL, intent(out):: ps(:, :) ! (iim + 1, jjm + 1) in Pa

    ! Local variables: 
    INTEGER iq
    REAL, allocatable:: tab_cntrl(:) ! tableau des param\`etres du run
    INTEGER ierr, ncid, varid

    !-----------------------------------------------------------------------

    print *, "Call sequence information: dynetat0"

    call assert((/size(ucov, 1), size(vcov, 1), size(masse, 1), size(ps, 1), &
         size(q, 1), size(teta, 1)/) == iim + 1, "dynetat0 iim")
    call assert((/size(ucov, 2), size(vcov, 2) + 1, size(masse, 2), &
         size(ps, 2), size(q, 2), size(teta, 2)/) == jjm + 1, "dynetat0 jjm")
    call assert((/size(vcov, 3), size(ucov, 3), size(teta, 3), size(q, 3), &
         size(masse, 3)/) == llm, "dynetat0 llm")
    call assert(size(q, 4) == nqmx, "dynetat0 q nqmx")

    allocate(rlatu(jjm + 1))
    allocate(rlatv(jjm))
    allocate(rlonu(iim + 1))
    allocate(rlonv(iim + 1))
    allocate(xprimu(iim + 1), xprimv(iim + 1))
    allocate(xprimm025(iim + 1), xprimp025(iim + 1))
    allocate(rlatu1(jjm), rlatu2(jjm), yprimu1(jjm), yprimu2(jjm))

    ! Fichier \'etat initial :
    call nf95_open("start.nc", NF90_NOWRITE, ncid)

    call nf95_inq_varid(ncid, "controle", varid)
    call NF95_Gw_VAR(ncid, varid, tab_cntrl)

    etot0 = tab_cntrl(13)
    ptot0 = tab_cntrl(14)
    ztot0 = tab_cntrl(15)
    stot0 = tab_cntrl(16)
    ang0 = tab_cntrl(17)

    if (raz_date) then
       print *, 'Resetting the date.'
       day_ini = day_ref
       itau_dyn = 0
    else
       itau_dyn = tab_cntrl(31)
       day_ini = tab_cntrl(30)
    end if

    print *, "day_ini = ", day_ini

    call NF95_INQ_VARID (ncid, "rlonu", varid)
    call NF95_GET_VAR(ncid, varid, rlonu)

    call NF95_INQ_VARID (ncid, "rlatu", varid)
    call NF95_GET_VAR(ncid, varid, rlatu)

    call NF95_INQ_VARID (ncid, "rlonv", varid)
    call NF95_GET_VAR(ncid, varid, rlonv)

    call NF95_INQ_VARID (ncid, "rlatv", varid)
    call NF95_GET_VAR(ncid, varid, rlatv)

    CALL nf95_inq_varid(ncid, 'xprimu', varid)
    CALL nf95_get_var(ncid, varid, xprimu)

    CALL nf95_inq_varid(ncid, 'xprimv', varid)
    CALL nf95_get_var(ncid, varid, xprimv)

    CALL nf95_inq_varid(ncid, 'xprimm025', varid)
    CALL nf95_get_var(ncid, varid, xprimm025)

    CALL nf95_inq_varid(ncid, 'xprimp025', varid)
    CALL nf95_get_var(ncid, varid, xprimp025)

    call NF95_INQ_VARID (ncid, "rlatu1", varid)
    call NF95_GET_VAR(ncid, varid, rlatu1)

    call NF95_INQ_VARID (ncid, "rlatu2", varid)
    call NF95_GET_VAR(ncid, varid, rlatu2)

    CALL nf95_inq_varid(ncid, 'yprimu1', varid)
    CALL nf95_get_var(ncid, varid, yprimu1)

    CALL nf95_inq_varid(ncid, 'yprimu2', varid)
    CALL nf95_get_var(ncid, varid, yprimu2)

    call NF95_INQ_VARID (ncid, "phis", varid)
    call NF95_GET_VAR(ncid, varid, phis)

    call NF95_INQ_VARID (ncid, "ucov", varid)
    call NF95_GET_VAR(ncid, varid, ucov)

    call NF95_INQ_VARID (ncid, "vcov", varid)
    call NF95_GET_VAR(ncid, varid, vcov)

    call NF95_INQ_VARID (ncid, "teta", varid)
    call NF95_GET_VAR(ncid, varid, teta)

    DO iq = 1, nqmx
       call NF95_INQ_VARID(ncid, tname(iq), varid, ierr)
       IF (ierr == NF90_NOERR) THEN
          call NF95_GET_VAR(ncid, varid, q(:, :, :, iq))
       ELSE
          PRINT *, 'dynetat0: "' // tname(iq) // '" not found, ' // &
               "setting it to zero..."
          q(:, :, :, iq) = 0.
       ENDIF
    ENDDO

    call NF95_INQ_VARID (ncid, "masse", varid)
    call NF95_GET_VAR(ncid, varid, masse)

    call NF95_INQ_VARID (ncid, "ps", varid)
    call NF95_GET_VAR(ncid, varid, ps)
    ! Check that there is a single value at each pole:
    call assert(ps(1, 1) == ps(2:, 1), "dynetat0 ps north pole")
    call assert(ps(1, jjm + 1) == ps(2:, jjm + 1), "dynetat0 ps south pole")

    call NF95_CLOSE(ncid)

  END SUBROUTINE dynetat0

  !********************************************************************

  SUBROUTINE fyhyp

    ! From LMDZ4/libf/dyn3d/fyhyp.F, version 1.2, 2005/06/03 09:11:32

    ! Author: P. Le Van, from analysis by R. Sadourny 

    ! Define rlatu, rlatv, rlatu2, yprimu2, rlatu1, yprimu1.
    
    ! Calcule les latitudes et dérivées dans la grille du GCM pour une
    ! fonction f(y) à dérivée tangente hyperbolique.

    ! Il vaut mieux avoir : grossismy * dzoom < pi / 2

    use jumble, only: rad_to_deg

    use coefpoly_m, only: coefpoly, a0, a1, a2, a3
    USE dimensions, only: jjm
    use dynetat0_chosen_m, only: clat, grossismy, dzoomy, tauy
    use heavyside_m, only: heavyside

    ! Local: 

    INTEGER, PARAMETER:: nmax2 = 2 * nmax
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
    DOUBLE PRECISION ypn
    DOUBLE PRECISION, save::deply, y00

    INTEGER i, j, it, ik, iter, jlat, jjpn
    INTEGER, save:: jpn
    DOUBLE PRECISION yi2, heavyy0, heavyy0m
    DOUBLE PRECISION fa(0:nmax2), fb(0:nmax2)
    REAL y0min, y0max

    !-------------------------------------------------------------------

    print *, "Call sequence information: fyhyp"

    allocate(rlatu(jjm + 1))
    allocate(rlatv(jjm))
    allocate(rlatu1(jjm), rlatu2(jjm), yprimu1(jjm), yprimu2(jjm))

    pi = 2.*asin(1.)
    pis2 = pi/2.
    pisjm = pi/real(jjm)
    epsilon = 1e-3
    dzoom = dzoomy*pi

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
               yt(it), yt(it + 1))

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

    print *, "Minimum latitude step:", minval(ylat(:jjm)) * rad_to_deg, &
         "degrees"
    print *, "Maximum latitude step:", maxval(ylat(:jjm)) * rad_to_deg, &
         "degrees"

    rlatu(1) = pi / 2.
    rlatu(jjm + 1) = -rlatu(1)

  END SUBROUTINE fyhyp

  !********************************************************************

  SUBROUTINE fxhyp

    ! From LMDZ4/libf/dyn3d/fxhyp.F, version 1.2, 2005/06/03 09:11:32
    ! Author: P. Le Van, from formulas by R. Sadourny

    ! Compute xprimm025, rlonv, xprimv, rlonu, xprimu, xprimp025.
    
    ! Calcule les longitudes et dérivées dans la grille du GCM pour
    ! une fonction $x_f(\tilde x)$ à dérivée tangente hyperbolique.

    ! Il vaut mieux avoir : grossismx $\times$ delta < pi

    ! Le premier point scalaire pour une grille regulière (grossismx =
    ! 1) avec clon = 0 est à - 180 degrés.

    ! Libraries:
    use jumble, only: pi, pi_d, twopi, twopi_d, arth, assert, rad_to_deg

    USE dimensions, ONLY: iim
    use dynetat0_chosen_m, only: clon, grossismx, dzoomx, taux
    use invert_zoom_x_m, only: invert_zoom_x
    use principal_cshift_m, only: principal_cshift
    use tanh_cautious_m, only: tanh_cautious

    ! Local:
    real rlonm025(iim + 1), rlonp025(iim + 1), d_rlonv(iim)
    REAL delta, h
    DOUBLE PRECISION, dimension(0:nmax):: xtild, fhyp, G, Xf, ffdx
    DOUBLE PRECISION beta
    INTEGER i, is2
    DOUBLE PRECISION xmoy(nmax), fxm(nmax)

    !----------------------------------------------------------------------

    print *, "Call sequence information: fxhyp"

    allocate(rlonu(iim + 1))
    allocate(rlonv(iim + 1))
    allocate(xprimu(iim + 1), xprimv(iim + 1))
    allocate(xprimm025(iim + 1), xprimp025(iim + 1))

    if (grossismx == 1.) then
       h = twopi / iim

       xprimm025(:iim) = h
       xprimp025(:iim) = h
       xprimv(:iim) = h
       xprimu(:iim) = h

       rlonv(:iim) = arth(- pi + clon, h, iim)
       rlonm025(:iim) = rlonv(:iim) - 0.25 * h
       rlonp025(:iim) = rlonv(:iim) + 0.25 * h
       rlonu(:iim) = rlonv(:iim) + 0.5 * h
    else
       delta = dzoomx * twopi_d
       xtild = arth(0d0, pi_d / nmax, nmax + 1)
       forall (i = 1:nmax) xmoy(i) = 0.5d0 * (xtild(i-1) + xtild(i))

       ! Compute fhyp:
       fhyp(1:nmax - 1) = tanh_cautious(taux * (delta / 2d0 &
            - xtild(1:nmax - 1)), xtild(1:nmax - 1) &
            * (pi_d - xtild(1:nmax - 1)))
       fhyp(0) = 1d0
       fhyp(nmax) = -1d0

       fxm = tanh_cautious(taux * (delta / 2d0 - xmoy), xmoy * (pi_d - xmoy))

       ! Compute \int_0 ^{\tilde x} F:

       ffdx(0) = 0d0

       DO i = 1, nmax
          ffdx(i) = ffdx(i - 1) + fxm(i) * (xtild(i) - xtild(i-1))
       END DO

       print *, "ffdx(nmax) = ", ffdx(nmax)
       beta = (pi_d - grossismx * ffdx(nmax)) / (pi_d - ffdx(nmax))
       print *, "beta = ", beta

       IF (2d0 * beta - grossismx <= 0d0) THEN
          print *, 'Bad choice of grossismx, taux, dzoomx.'
          print *, 'Decrease dzoomx or grossismx.'
          STOP 1
       END IF

       G = beta + (grossismx - beta) * fhyp

       Xf(:nmax - 1) = beta * xtild(:nmax - 1) + (grossismx - beta) &
            * ffdx(:nmax - 1)
       Xf(nmax) = pi_d

       call invert_zoom_x(beta, xf, xtild, G, rlonm025(:iim), xprimm025(:iim), &
            xuv = - 0.25d0)
       call invert_zoom_x(beta, xf, xtild, G, rlonv(:iim), xprimv(:iim), &
            xuv = 0d0)
       call invert_zoom_x(beta, xf, xtild, G, rlonu(:iim), xprimu(:iim), &
            xuv = 0.5d0)
       call invert_zoom_x(beta, xf, xtild, G, rlonp025(:iim), xprimp025(:iim), &
            xuv = 0.25d0)
    end if

    is2 = 0

    IF (MINval(rlonm025(:iim)) < - pi - 0.1 &
         .or. MAXval(rlonm025(:iim)) > pi + 0.1) THEN
       IF (clon <= 0.) THEN
          is2 = 1

          do while (rlonm025(is2) < - pi .and. is2 < iim)
             is2 = is2 + 1
          end do

          call assert(rlonm025(is2) >= - pi, &
               "fxhyp -- rlonm025 should be >= - pi")
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
    print *, "Minimum longitude step:", MINval(d_rlonv) * rad_to_deg, "degrees"
    print *, "Maximum longitude step:", MAXval(d_rlonv) * rad_to_deg, "degrees"

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

end module dynetat0_m
