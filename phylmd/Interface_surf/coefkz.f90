module coefkz_m

  IMPLICIT none

contains

  SUBROUTINE coefkz(nsrf, paprs, pplay, ts, u, v, t, q, zgeop, coefm, coefh)

    ! Authors: F. Hourdin, M. Forichon, Z. X. Li (LMD/CNRS)
    ! Date: September 22nd, 1993

    ! Objet : calculer les coefficients d'échange turbulent dans
    ! l'atmosphère.

    USE clesphys, ONLY: ksta, ksta_ter
    USE conf_phys_m, ONLY: iflag_pbl
    USE dimphy, ONLY: klev
    USE fcttre, ONLY: foede, foeew
    USE indicesol, ONLY: is_oce
    USE suphec_m, ONLY: rcpd, rd, retv, rg, rkappa, rlstt, rlvtt, rtt
    USE yoethf_m, ONLY: r2es, r5ies, r5les, rvtmp2

    integer, intent(in):: nsrf ! indicateur de la nature du sol

    REAL, intent(in):: paprs(:, :) ! (knon, klev + 1)
    ! pression a chaque intercouche (en Pa)

    real, intent(in):: pplay(:, :) ! (knon, klev)
    ! pression au milieu de chaque couche (en Pa)

    REAL, intent(in):: ts(:) ! (knon) temperature du sol (en Kelvin)
    REAL, intent(in):: u(:, :), v(:, :) ! (knon, klev) wind
    REAL, intent(in):: t(:, :) ! (knon, klev) temperature (K)
    real, intent(in):: q(:, :) ! (knon, klev) vapeur d'eau (kg / kg)
    REAL, intent(in):: zgeop(:, :) ! (knon, klev)
    REAL, intent(out):: coefm(:, 2:) ! (knon, 2:klev) coefficient, vitesse

    real, intent(out):: coefh(:, 2:) ! (knon, 2:klev)
    ! coefficient, chaleur et humidité

    ! Local:

    INTEGER knon ! nombre de points a traiter

    INTEGER itop(size(ts)) ! (knon)
    ! numero de couche du sommet de la couche limite

    ! Quelques constantes et options:

    REAL, PARAMETER:: cepdu2 =0.1**2
    REAL, PARAMETER:: ratqs = 0.05 ! largeur de distribution de vapeur d'eau
    REAL, PARAMETER:: ric = 0.4 ! nombre de Richardson critique
    REAL, PARAMETER:: prandtl = 0.4

    REAL kstable ! diffusion minimale (situation stable)
    REAL, PARAMETER:: mixlen = 35. ! constante contrôlant longueur de mélange
    INTEGER i, k
    REAL zmgeom(size(ts))
    REAL ri(size(ts))
    REAL l2(size(ts))
    REAL zdphi, zdu2, ztvd, ztvu, cdn
    REAL zt, zq, zcvm5, zcor, zqs, zfr, zdqs
    logical zdelta
    REAL gamt(2:klev) ! contre-gradient pour la chaleur sensible: Kelvin / metre

    !--------------------------------------------------------------------

    knon = size(ts)

    ! Prescrire la valeur de contre-gradient
    if (iflag_pbl == 1) then
       DO k = 3, klev
          gamt(k) = - 1E-3
       ENDDO
       gamt(2) = - 2.5E-3
    else
       DO k = 2, klev
          gamt(k) = 0.
       ENDDO
    ENDIF

    kstable = merge(ksta, ksta_ter, nsrf == is_oce)

    ! Calculer les coefficients turbulents dans l'atmosphere

    itop = klev

    DO k = 2, klev
       DO i = 1, knon
          zdu2 = MAX(cepdu2, (u(i, k) - u(i, k - 1))**2 &
               + (v(i, k) - v(i, k - 1))**2)
          zmgeom(i) = zgeop(i, k) - zgeop(i, k - 1)
          zdphi = zmgeom(i) / 2.
          zt = (t(i, k) + t(i, k - 1)) * 0.5
          zq = (q(i, k) + q(i, k - 1)) * 0.5

          ! calculer Qs et dQs / dT:
          zdelta = RTT >=zt
          zcvm5 = merge(R5IES * RLSTT, R5LES * RLVTT, zdelta) / RCPD &
               / (1. + RVTMP2 * zq)
          zqs = R2ES * FOEEW(zt, zdelta) / pplay(i, k)
          zqs = MIN(0.5, zqs)
          zcor = 1. / (1. - RETV * zqs)
          zqs = zqs * zcor
          zdqs = FOEDE(zt, zdelta, zcvm5, zqs, zcor)

          ! calculer la fraction nuageuse (processus humide):
          zfr = (zq + ratqs * zq - zqs) / (2. * ratqs * zq)
          zfr = MAX(0., MIN(1., zfr))

          ! calculer le nombre de Richardson:
          ztvd = (t(i, k) + zdphi / RCPD / (1. + RVTMP2 * zq) * ((1. - zfr) &
               + zfr * (1. + RLVTT * zqs / RD / zt) / (1. + zdqs))) * (1. &
               + RETV * q(i, k))
          ztvu = (t(i, k - 1) - zdphi / RCPD / (1. + RVTMP2 * zq) * ((1. &
               - zfr) + zfr * (1. + RLVTT * zqs / RD / zt) / (1. + zdqs))) &
               * (1. + RETV * q(i, k - 1))
          ri(i) = zmgeom(i) * (ztvd - ztvu) / (zdu2 * 0.5 * (ztvd + ztvu))
          ri(i) = ri(i) &
          cdn = SQRT(zdu2) / zmgeom(i) * RG

          l2(i) = (mixlen * MAX(0.0, (paprs(i, k) - paprs(i, itop(i) + 1)) &
               /(paprs(i, 2) - paprs(i, itop(i) + 1))))**2
          coefm(i, k) = sqrt(max(cdn**2 * (ric - ri(i)) / ric, kstable))
          coefm(i, k) = l2(i) * coefm(i, k)
               + zmgeom(i) * zmgeom(i) / RG * gamt(k) &
               * (paprs(i, k) / 101325.)**RKAPPA &
               / (zdu2 * 0.5 * (ztvd + ztvu))

          ! Finalement, les coefficients d'\'echange sont obtenus:
          coefh(i, k) = coefm(i, k) / prandtl ! h et m different
       ENDDO
    ENDDO

    ! Au-delà du sommet, pas de diffusion turbulente :
    forall (i = 1: knon)
       coefh(i, itop(i) + 1:) = 0.
       coefm(i, itop(i) + 1:) = 0.
    END forall

  END SUBROUTINE coefkz

end module coefkz_m
