module clqh_m

  IMPLICIT none

contains

  SUBROUTINE clqh(dtime, julien, debut, nisurf, knindex, tsoil, qsol, rmu0, &
       rugos, rugoro, u1lay, v1lay, coef, tq_cdrag, t, q, ts, paprs, pplay, &
       delp, radsol, albedo, snow, qsurf, precip_rain, precip_snow, fluxlat, &
       pctsrf_new_sic, agesno, d_t, d_q, d_ts, z0_new, flux_t, flux_q, &
       dflux_s, dflux_l, fqcalving, ffonte, run_off_lic_0)

    ! Author: Z. X. Li (LMD/CNRS)
    ! Date: 1993 Aug. 18th
    ! Objet : diffusion verticale de "q" et de "h"

    USE conf_phys_m, ONLY: iflag_pbl
    USE dimphy, ONLY: klev, klon
    USE interfsurf_hq_m, ONLY: interfsurf_hq
    USE suphec_m, ONLY: rcpd, rd, rg, rkappa

    REAL, intent(in):: dtime ! intervalle du temps (s)
    integer, intent(in):: julien ! jour de l'annee en cours
    logical, intent(in):: debut
    integer, intent(in):: nisurf
    integer, intent(in):: knindex(:) ! (knon)
    REAL, intent(inout):: tsoil(:, :) ! (knon, nsoilmx)

    REAL, intent(inout):: qsol(:) ! (knon)
    ! column-density of water in soil, in kg m-2

    real, intent(in):: rmu0(klon) ! cosinus de l'angle solaire zenithal
    real, intent(in):: rugos(:) ! (knon) rugosite
    REAL, intent(in):: rugoro(:) ! (knon)

    REAL, intent(in):: u1lay(:), v1lay(:) ! (knon)
    ! vitesse de la 1ere couche (m / s)

    REAL, intent(in):: coef(:, 2:) ! (knon, 2:klev)
    ! Le coefficient d'echange (m**2 / s) multiplie par le cisaillement
    ! du vent (dV / dz)

    REAL, intent(in):: tq_cdrag(:) ! (knon) sans unite

    REAL, intent(in):: t(:, :) ! (knon, klev) temperature (K)
    REAL, intent(in):: q(:, :) ! (knon, klev) humidite specifique (kg / kg)
    REAL, intent(in):: ts(:) ! (knon) temperature du sol (K)

    REAL, intent(in):: paprs(:, :) ! (knon, klev + 1)
    ! pression a inter-couche (Pa)

    REAL, intent(in):: pplay(:, :) ! (knon, klev)
    ! pression au milieu de couche (Pa)

    REAL, intent(in):: delp(:, :) ! (knon, klev)
    ! epaisseur de couche en pression (Pa)

    REAL, intent(in):: radsol(:) ! (knon)
    ! rayonnement net au sol (Solaire + IR) W / m2

    REAL, intent(inout):: albedo(:) ! (knon) albedo de la surface
    REAL, intent(inout):: snow(:) ! (knon) ! hauteur de neige

    REAL, intent(out):: qsurf(:) ! (knon)
    ! humidite de l'air au dessus de la surface

    real, intent(in):: precip_rain(klon)
    ! liquid water mass flux (kg / m2 / s), positive down

    real, intent(in):: precip_snow(klon)
    ! solid water mass flux (kg / m2 / s), positive down

    real, intent(out):: fluxlat(:) ! (knon)
    real, intent(in):: pctsrf_new_sic(:) ! (klon)
    REAL, intent(inout):: agesno(:) ! (knon)
    REAL, intent(out):: d_t(:, :) ! (knon, klev) incrementation de "t"
    REAL, intent(out):: d_q(:, :) ! (knon, klev) incrementation de "q"
    REAL, intent(out):: d_ts(:) ! (knon) variation of surface temperature
    real, intent(out):: z0_new(:) ! (knon)

    REAL, intent(out):: flux_t(:) ! (knon)
    ! (diagnostic) flux de chaleur sensible (Cp T) à la surface,
    ! positif vers le bas, W / m2

    REAL, intent(out):: flux_q(:) ! (knon)
    ! flux de la vapeur d'eau à la surface, en kg / (m**2 s)

    REAL, intent(out):: dflux_s(:) ! (knon) derivee du flux sensible dF / dTs
    REAL, intent(out):: dflux_l(:) ! (knon) derivee du flux latent dF / dTs

    REAL, intent(out):: fqcalving(:) ! (knon)
    ! Flux d'eau "perdue" par la surface et n\'ecessaire pour que limiter la
    ! hauteur de neige, en kg / m2 / s

    REAL ffonte(klon)
    ! Flux thermique utiliser pour fondre la neige

    REAL run_off_lic_0(klon)! runof glacier au pas de temps precedent

    ! Local:
    INTEGER knon
    REAL evap(size(knindex)) ! (knon) evaporation au sol
    INTEGER i, k
    REAL, dimension(size(knindex), klev):: cq, dq, ch, dh ! (knon, klev)
    REAL buf1(size(knindex)), buf2(size(knindex))
    REAL zx_coef(size(knindex), 2:klev) ! (knon, 2:klev)
    REAL h(size(knindex), klev) ! (knon, klev) enthalpie potentielle
    REAL local_q(size(knindex), klev) ! (knon, klev)

    REAL psref(size(knindex)) ! (knon)
    ! pression de reference pour temperature potentielle

    REAL pkf(size(knindex), klev) ! (knon, klev)

    REAL gamt(size(knindex), 2:klev) ! (knon, 2:klev)
    ! contre-gradient pour la chaleur sensible, en K m-1

    REAL gamah(size(knindex), 2:klev) ! (knon, 2:klev)
    real tsurf_new(size(knindex)) ! (knon)

    !----------------------------------------------------------------

    knon = size(knindex)

    if (iflag_pbl == 1) then
       gamt(:, 2) = - 2.5e-3
       gamt(:, 3:)= - 1e-3
    else
       gamt = 0.
    endif

    psref = paprs(:, 1) ! pression de reference est celle au sol
    forall (k = 1:klev) pkf(:, k) = (psref / pplay(:, k))**RKAPPA
    h = RCPD * t * pkf

    ! Convertir les coefficients en variables convenables au calcul:
    forall (k = 2:klev) zx_coef(:, k) = coef(:, k) * RG &
         / (pplay(:, k - 1) - pplay(:, k)) &
         * (paprs(:, k) * 2 / (t(:, k) + t(:, k - 1)) / RD)**2 * dtime * RG

    ! Preparer les flux lies aux contre-gardients
    forall (k = 2:klev) gamah(:, k) = gamt(:, k) * (RD * (t(:, k - 1) &
         + t(:, k)) / 2. / RG / paprs(:, k) * (pplay(:, k - 1) - pplay(:, k))) &
         * RCPD * (psref / paprs(:, k))**RKAPPA

    DO i = 1, knon
       buf1(i) = zx_coef(i, klev) + delp(i, klev)
       cq(i, klev) = q(i, klev) * delp(i, klev) / buf1(i)
       dq(i, klev) = zx_coef(i, klev) / buf1(i)

       buf2(i) = delp(i, klev) / pkf(i, klev) + zx_coef(i, klev)
       ch(i, klev) = (h(i, klev) / pkf(i, klev) * delp(i, klev) &
            - zx_coef(i, klev) * gamah(i, klev)) / buf2(i)
       dh(i, klev) = zx_coef(i, klev) / buf2(i)
    ENDDO

    DO k = klev - 1, 2, - 1
       DO i = 1, knon
          buf1(i) = delp(i, k) + zx_coef(i, k) &
               + zx_coef(i, k + 1) * (1. - dq(i, k + 1))
          cq(i, k) = (q(i, k) * delp(i, k) &
               + zx_coef(i, k + 1) * cq(i, k + 1)) / buf1(i)
          dq(i, k) = zx_coef(i, k) / buf1(i)

          buf2(i) = delp(i, k) / pkf(i, k) + zx_coef(i, k) &
               + zx_coef(i, k + 1) * (1. - dh(i, k + 1))
          ch(i, k) = (h(i, k) / pkf(i, k) * delp(i, k) &
               + zx_coef(i, k + 1) * ch(i, k + 1) &
               + zx_coef(i, k + 1) * gamah(i, k + 1) &
               - zx_coef(i, k) * gamah(i, k)) / buf2(i)
          dh(i, k) = zx_coef(i, k) / buf2(i)
       ENDDO
    ENDDO

    DO i = 1, knon
       buf1(i) = delp(i, 1) + zx_coef(i, 2) * (1. - dq(i, 2))
       cq(i, 1) = (q(i, 1) * delp(i, 1) + zx_coef(i, 2) * cq(i, 2)) / buf1(i)
       dq(i, 1) = - 1. * RG / buf1(i)

       buf2(i) = delp(i, 1) / pkf(i, 1) + zx_coef(i, 2) * (1. - dh(i, 2))
       ch(i, 1) = (h(i, 1) / pkf(i, 1) * delp(i, 1) &
            + zx_coef(i, 2) * (gamah(i, 2) + ch(i, 2))) / buf2(i)
       dh(i, 1) = - 1. * RG / buf2(i)
    ENDDO

    CALL interfsurf_hq(dtime, julien, rmu0, nisurf, knindex, debut, tsoil, &
         qsol, u1lay, v1lay, t(:, 1), q(:, 1), tq_cdrag(:knon), ch(:, 1), &
         cq(:, 1), dh(:, 1), dq(:, 1), precip_rain, precip_snow, rugos, &
         rugoro, snow, qsurf, ts, pplay(:, 1), psref, radsol, evap, flux_t, &
         fluxlat, dflux_l, dflux_s, tsurf_new, albedo, z0_new, pctsrf_new_sic, &
         agesno, fqcalving, ffonte, run_off_lic_0)

    flux_q = - evap
    d_ts = tsurf_new - ts

    h(:, 1) = ch(:, 1) + dh(:, 1) * flux_t * dtime
    local_q(:, 1) = cq(:, 1) + dq(:, 1) * flux_q * dtime

    DO k = 2, klev
       h(:, k) = ch(:, k) + dh(:, k) * h(:, k - 1)
       local_q(:, k) = cq(:, k) + dq(:, k) * local_q(:, k - 1)
    ENDDO

    d_t = h / pkf / RCPD - t
    d_q = local_q - q

  END SUBROUTINE clqh

end module clqh_m
