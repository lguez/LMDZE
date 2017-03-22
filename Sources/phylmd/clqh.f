module clqh_m

  IMPLICIT none

contains

  SUBROUTINE clqh(dtime, jour, debut, nisurf, knindex, tsoil, qsol, rmu0, &
       rugos, rugoro, u1lay, v1lay, coef, t, q, ts, paprs, pplay, delp, &
       radsol, albedo, snow, qsurf, precip_rain, precip_snow, fder, fluxlat, &
       pctsrf_new_sic, agesno, d_t, d_q, d_ts, z0_new, flux_t, flux_q, &
       dflux_s, dflux_l, fqcalving, ffonte, run_off_lic_0)

    ! Author: Z. X. Li (LMD/CNRS)
    ! Date: 1993/08/18
    ! Objet : diffusion verticale de "q" et de "h"

    USE conf_phys_m, ONLY: iflag_pbl
    USE dimphy, ONLY: klev, klon
    USE interfsurf_hq_m, ONLY: interfsurf_hq
    USE suphec_m, ONLY: rcpd, rd, rg, rkappa

    REAL, intent(in):: dtime ! intervalle du temps (s)
    integer, intent(in):: jour ! jour de l'annee en cours
    logical, intent(in):: debut
    integer, intent(in):: nisurf
    integer, intent(in):: knindex(:) ! (knon)
    REAL, intent(inout):: tsoil(:, :) ! (knon, nsoilmx)

    REAL, intent(inout):: qsol(klon)
    ! column-density of water in soil, in kg m-2

    real, intent(in):: rmu0(klon) ! cosinus de l'angle solaire zenithal
    real rugos(klon) ! rugosite
    REAL rugoro(klon)
    REAL u1lay(klon) ! vitesse u de la 1ere couche (m / s)
    REAL v1lay(klon) ! vitesse v de la 1ere couche (m / s)

    REAL, intent(in):: coef(:, :) ! (knon, klev)
    ! Le coefficient d'echange (m**2 / s) multiplie par le cisaillement
    ! du vent (dV / dz). La premiere valeur indique la valeur de Cdrag
    ! (sans unite).

    REAL t(klon, klev) ! temperature (K)
    REAL q(klon, klev) ! humidite specifique (kg / kg)
    REAL, intent(in):: ts(:) ! (knon) temperature du sol (K)
    REAL paprs(klon, klev + 1) ! pression a inter-couche (Pa)
    REAL pplay(klon, klev) ! pression au milieu de couche (Pa)
    REAL delp(klon, klev) ! epaisseur de couche en pression (Pa)
    REAL radsol(klon) ! ray. net au sol (Solaire + IR) W / m2
    REAL, intent(inout):: albedo(:) ! (knon) albedo de la surface
    REAL, intent(inout):: snow(klon) ! hauteur de neige
    REAL qsurf(klon) ! humidite de l'air au dessus de la surface

    real, intent(in):: precip_rain(klon)
    ! liquid water mass flux (kg / m2 / s), positive down

    real, intent(in):: precip_snow(klon)
    ! solid water mass flux (kg / m2 / s), positive down

    real, intent(inout):: fder(klon)
    real, intent(out):: fluxlat(:) ! (knon)
    real, intent(in):: pctsrf_new_sic(:) ! (klon)
    REAL, intent(inout):: agesno(:) ! (knon)
    REAL d_t(klon, klev) ! incrementation de "t"
    REAL d_q(klon, klev) ! incrementation de "q"
    REAL, intent(out):: d_ts(:) ! (knon) incr\'ementation de "ts"
    real z0_new(klon)

    REAL, intent(out):: flux_t(:) ! (knon)
    ! (diagnostic) flux de chaleur sensible (Cp T) à la surface,
    ! positif vers le bas, W / m2

    REAL, intent(out):: flux_q(:) ! (knon)
    ! flux de la vapeur d'eau à la surface, en kg / (m**2 s)

    REAL dflux_s(klon) ! derivee du flux sensible dF / dTs
    REAL dflux_l(klon) ! derivee du flux latent dF / dTs

    ! Flux d'eau "perdue" par la surface et n\'ecessaire pour que limiter la
    ! hauteur de neige, en kg / m2 / s
    REAL fqcalving(klon)

    ! Flux thermique utiliser pour fondre la neige
    REAL ffonte(klon)

    REAL run_off_lic_0(klon)! runof glacier au pas de temps precedent

    ! Local:

    INTEGER knon
    REAL evap(size(knindex)) ! (knon) evaporation au sol

    INTEGER i, k
    REAL zx_cq(klon, klev)
    REAL zx_dq(klon, klev)
    REAL zx_ch(klon, klev)
    REAL zx_dh(klon, klev)
    REAL zx_buf1(klon)
    REAL zx_buf2(klon)
    REAL zx_coef(klon, klev)
    REAL local_h(klon, klev) ! enthalpie potentielle
    REAL local_q(klon, klev)
    REAL psref(klon) ! pression de reference pour temperature potent.
    REAL zx_pkh(klon, klev), zx_pkf(klon, klev)

    ! contre-gradient pour la vapeur d'eau: (kg / kg) / metre
    REAL gamq(klon, 2:klev)
    ! contre-gradient pour la chaleur sensible: Kelvin / metre
    REAL gamt(klon, 2:klev)
    REAL z_gamaq(klon, 2:klev), z_gamah(klon, 2:klev)
    REAL zdelz

    real temp_air(klon), spechum(klon)
    real tq_cdrag(klon), petAcoef(klon), peqAcoef(klon)
    real petBcoef(klon), peqBcoef(klon)
    real p1lay(klon)

    real tsurf_new(size(knindex)) ! (knon)
    real zzpk

    !----------------------------------------------------------------

    knon = size(knindex)

    if (iflag_pbl == 1) then
       do k = 3, klev
          do i = 1, knon
             gamq(i, k)= 0.0
             gamt(i, k)= - 1.0e-03
          enddo
       enddo
       do i = 1, knon
          gamq(i, 2) = 0.0
          gamt(i, 2) = - 2.5e-03
       enddo
    else
       do k = 2, klev
          do i = 1, knon
             gamq(i, k) = 0.0
             gamt(i, k) = 0.0
          enddo
       enddo
    endif

    DO i = 1, knon
       psref(i) = paprs(i, 1) !pression de reference est celle au sol
    ENDDO
    DO k = 1, klev
       DO i = 1, knon
          zx_pkh(i, k) = (psref(i) / paprs(i, k))**RKAPPA
          zx_pkf(i, k) = (psref(i) / pplay(i, k))**RKAPPA
          local_h(i, k) = RCPD * t(i, k) * zx_pkf(i, k)
          local_q(i, k) = q(i, k)
       ENDDO
    ENDDO

    ! Convertir les coefficients en variables convenables au calcul:

    DO k = 2, klev
       DO i = 1, knon
          zx_coef(i, k) = coef(i, k) * RG / (pplay(i, k - 1) - pplay(i, k)) &
               * (paprs(i, k) * 2 / (t(i, k) + t(i, k - 1)) / RD)**2
          zx_coef(i, k) = zx_coef(i, k) * dtime * RG
       ENDDO
    ENDDO

    ! Preparer les flux lies aux contre-gardients

    DO k = 2, klev
       DO i = 1, knon
          zdelz = RD * (t(i, k - 1) + t(i, k)) / 2.0 / RG / paprs(i, k) &
               * (pplay(i, k - 1) - pplay(i, k))
          z_gamaq(i, k) = gamq(i, k) * zdelz
          z_gamah(i, k) = gamt(i, k) * zdelz * RCPD * zx_pkh(i, k)
       ENDDO
    ENDDO
    DO i = 1, knon
       zx_buf1(i) = zx_coef(i, klev) + delp(i, klev)
       zx_cq(i, klev) = (local_q(i, klev) * delp(i, klev) &
            - zx_coef(i, klev) * z_gamaq(i, klev)) / zx_buf1(i)
       zx_dq(i, klev) = zx_coef(i, klev) / zx_buf1(i)

       zzpk=(pplay(i, klev) / psref(i))**RKAPPA
       zx_buf2(i) = zzpk * delp(i, klev) + zx_coef(i, klev)
       zx_ch(i, klev) = (local_h(i, klev) * zzpk * delp(i, klev) &
            - zx_coef(i, klev) * z_gamah(i, klev)) / zx_buf2(i)
       zx_dh(i, klev) = zx_coef(i, klev) / zx_buf2(i)
    ENDDO
    DO k = klev - 1, 2, - 1
       DO i = 1, knon
          zx_buf1(i) = delp(i, k) + zx_coef(i, k) &
               + zx_coef(i, k + 1) * (1. - zx_dq(i, k + 1))
          zx_cq(i, k) = (local_q(i, k) * delp(i, k) &
               + zx_coef(i, k + 1) * zx_cq(i, k + 1) &
               + zx_coef(i, k + 1) * z_gamaq(i, k + 1) &
               - zx_coef(i, k) * z_gamaq(i, k)) / zx_buf1(i)
          zx_dq(i, k) = zx_coef(i, k) / zx_buf1(i)

          zzpk=(pplay(i, k) / psref(i))**RKAPPA
          zx_buf2(i) = zzpk * delp(i, k) + zx_coef(i, k) &
               + zx_coef(i, k + 1) * (1. - zx_dh(i, k + 1))
          zx_ch(i, k) = (local_h(i, k) * zzpk * delp(i, k) &
               + zx_coef(i, k + 1) * zx_ch(i, k + 1) &
               + zx_coef(i, k + 1) * z_gamah(i, k + 1) &
               - zx_coef(i, k) * z_gamah(i, k)) / zx_buf2(i)
          zx_dh(i, k) = zx_coef(i, k) / zx_buf2(i)
       ENDDO
    ENDDO

    DO i = 1, knon
       zx_buf1(i) = delp(i, 1) + zx_coef(i, 2) * (1. - zx_dq(i, 2))
       zx_cq(i, 1) = (local_q(i, 1) * delp(i, 1) &
            + zx_coef(i, 2) * (z_gamaq(i, 2) + zx_cq(i, 2))) / zx_buf1(i)
       zx_dq(i, 1) = - 1. * RG / zx_buf1(i)

       zzpk=(pplay(i, 1) / psref(i))**RKAPPA
       zx_buf2(i) = zzpk * delp(i, 1) + zx_coef(i, 2) * (1. - zx_dh(i, 2))
       zx_ch(i, 1) = (local_h(i, 1) * zzpk * delp(i, 1) &
            + zx_coef(i, 2) * (z_gamah(i, 2) + zx_ch(i, 2))) / zx_buf2(i)
       zx_dh(i, 1) = - 1. * RG / zx_buf2(i)
    ENDDO

    ! Appel \`a interfsurf (appel g\'en\'erique) routine d'interface
    ! avec la surface

    ! Initialisation
    petAcoef =0. 
    peqAcoef = 0.
    petBcoef =0.
    peqBcoef = 0.
    p1lay =0.

    petAcoef(1:knon) = zx_ch(1:knon, 1)
    peqAcoef(1:knon) = zx_cq(1:knon, 1)
    petBcoef(1:knon) = zx_dh(1:knon, 1)
    peqBcoef(1:knon) = zx_dq(1:knon, 1)
    tq_cdrag(1:knon) =coef(:knon, 1)
    temp_air(1:knon) =t(1:knon, 1)
    spechum(1:knon)=q(1:knon, 1)
    p1lay(1:knon) = pplay(1:knon, 1)

    CALL interfsurf_hq(dtime, jour, rmu0, nisurf, knon, knindex, debut, &
         tsoil, qsol, u1lay, v1lay, temp_air, spechum, tq_cdrag, petAcoef, &
         peqAcoef, petBcoef, peqBcoef, precip_rain, precip_snow, fder, rugos, &
         rugoro, snow, qsurf, ts, p1lay, psref, radsol, evap, flux_t, &
         fluxlat, dflux_l, dflux_s, tsurf_new, albedo, z0_new, &
         pctsrf_new_sic, agesno, fqcalving, ffonte, run_off_lic_0)

    flux_q = - evap
    d_ts = tsurf_new - ts

    ! Une fois qu'on a zx_h_ts, on peut faire l'it\'eration
    DO i = 1, knon
       local_h(i, 1) = zx_ch(i, 1) + zx_dh(i, 1) * flux_t(i) * dtime
       local_q(i, 1) = zx_cq(i, 1) + zx_dq(i, 1) * flux_q(i) * dtime
    ENDDO
    DO k = 2, klev
       DO i = 1, knon
          local_q(i, k) = zx_cq(i, k) + zx_dq(i, k) * local_q(i, k - 1)
          local_h(i, k) = zx_ch(i, k) + zx_dh(i, k) * local_h(i, k - 1)
       ENDDO
    ENDDO

    ! Calcul des tendances
    DO k = 1, klev
       DO i = 1, knon
          d_t(i, k) = local_h(i, k) / zx_pkf(i, k) / RCPD - t(i, k)
          d_q(i, k) = local_q(i, k) - q(i, k)
       ENDDO
    ENDDO

  END SUBROUTINE clqh

end module clqh_m
