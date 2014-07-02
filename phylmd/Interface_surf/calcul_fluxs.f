module calcul_fluxs_m

  implicit none

contains

  SUBROUTINE calcul_fluxs( klon, knon, nisurf, dtime,  &
       tsurf, p1lay, cal, beta, coef1lay, ps,  &
       precip_rain, precip_snow, snow, qsurf,  &
       radsol, dif_grnd, t1lay, q1lay, u1lay, v1lay,  &
       petAcoef, peqAcoef, petBcoef, peqBcoef,  &
       tsurf_new, evap, fluxlat, fluxsens, dflux_s, dflux_l)

    ! Cette routine calcule les fluxs en h et q à l'interface et une
    ! température de surface.

    ! L. Fairhead 4/2000

    ! input:
    ! knon nombre de points a traiter
    ! nisurf surface a traiter
    ! tsurf temperature de surface
    ! p1lay pression 1er niveau (milieu de couche)
    ! cal capacite calorifique du sol
    ! beta evap reelle
    ! coef1lay coefficient d'echange
    ! ps pression au sol
    ! precip_rain precipitations liquides
    ! precip_snow precipitations solides
    ! snow champs hauteur de neige
    ! runoff runoff en cas de trop plein
    ! petAcoef coeff. A de la resolution de la CL pour t
    ! peqAcoef coeff. A de la resolution de la CL pour q
    ! petBcoef coeff. B de la resolution de la CL pour t
    ! peqBcoef coeff. B de la resolution de la CL pour q
    ! radsol rayonnement net aus sol (LW + SW)
    ! dif_grnd coeff. diffusion vers le sol profond

    ! output:
    ! tsurf_new temperature au sol
    ! qsurf humidite de l'air au dessus du sol
    ! fluxsens flux de chaleur sensible
    ! fluxlat flux de chaleur latente
    ! dflux_s derivee du flux de chaleur sensible / Ts
    ! dflux_l derivee du flux de chaleur latente / Ts


    use indicesol
    use abort_gcm_m, only: abort_gcm
    use yoethf_m
    use fcttre, only: thermcep, foeew, qsats, qsatl, foede, dqsats, dqsatl
    use SUPHEC_M
    use interface_surf

    ! Parametres d'entree
    integer, intent(IN) :: knon, nisurf, klon
    real , intent(IN) :: dtime
    real, dimension(klon), intent(IN) :: petAcoef, peqAcoef
    real, dimension(klon), intent(IN) :: petBcoef, peqBcoef
    real, dimension(klon), intent(IN) :: ps, q1lay
    real, dimension(klon), intent(IN) :: tsurf, p1lay, cal, beta, coef1lay
    real, dimension(klon), intent(IN) :: precip_rain, precip_snow
    real, dimension(klon), intent(IN) :: radsol, dif_grnd
    real, dimension(klon), intent(IN) :: t1lay, u1lay, v1lay
    real, dimension(klon), intent(INOUT) :: snow, qsurf

    ! Parametres sorties
    real, dimension(klon), intent(OUT):: tsurf_new, evap, fluxsens, fluxlat
    real, dimension(klon), intent(OUT):: dflux_s, dflux_l

    ! Variables locales
    integer :: i
    real, dimension(klon) :: zx_mh, zx_nh, zx_oh
    real, dimension(klon) :: zx_mq, zx_nq, zx_oq
    real, dimension(klon) :: zx_pkh, zx_dq_s_dt, zx_qsat, zx_coef
    real, dimension(klon) :: zx_sl, zx_k1
    real, dimension(klon) :: zx_q_0 , d_ts
    real :: zdelta, zcvm5, zx_qs, zcor, zx_dq_s_dh
    real :: bilan_f, fq_fonte
    REAL :: subli, fsno
    REAL :: qsat_new, q1_new
    real, parameter :: t_grnd = 271.35, t_coup = 273.15
    !! PB temporaire en attendant mieux pour le modele de neige
    REAL, parameter :: chasno = 3.334E+05/(2.3867E+06*0.15)

    logical, save :: check = .false.
    character (len = 20) :: modname = 'calcul_fluxs'
    logical, save :: fonte_neige = .false.
    real, save :: max_eau_sol = 150.0
    character (len = 80) :: abort_message
    logical, save :: first = .true., second=.false.

    if (check) write(*, *)'Entree ', modname, ' surface = ', nisurf

    IF (check) THEN
       WRITE(*, *)' radsol (min, max)' &
            , MINVAL(radsol(1:knon)), MAXVAL(radsol(1:knon))
       !!CALL flush(6)
    ENDIF

    if (size(coastalflow) /= knon .AND. nisurf == is_ter) then
       write(*, *)'Bizarre, le nombre de points continentaux'
       write(*, *)'a change entre deux appels. J''arrete ...'
       abort_message='Pb run_off'
       call abort_gcm(modname, abort_message, 1)
    endif

    ! Traitement neige et humidite du sol

    ! Initialisation

    evap = 0.
    fluxsens=0.
    fluxlat=0.
    dflux_s = 0.
    dflux_l = 0.

    ! zx_qs = qsat en kg/kg

    DO i = 1, knon
       zx_pkh(i) = (ps(i)/ps(i))**RKAPPA
       IF (thermcep) THEN
          zdelta=MAX(0., SIGN(1., rtt-tsurf(i)))
          zcvm5 = R5LES*RLVTT*(1.-zdelta) + R5IES*RLSTT*zdelta
          zcvm5 = zcvm5 / RCPD / (1.0+RVTMP2*q1lay(i))
          zx_qs= r2es * FOEEW(tsurf(i), zdelta)/ps(i)
          zx_qs=MIN(0.5, zx_qs)
          zcor=1./(1.-retv*zx_qs)
          zx_qs=zx_qs*zcor
          zx_dq_s_dh = FOEDE(tsurf(i), zdelta, zcvm5, zx_qs, zcor) &
               /RLVTT / zx_pkh(i)
       ELSE
          IF (tsurf(i).LT.t_coup) THEN
             zx_qs = qsats(tsurf(i)) / ps(i)
             zx_dq_s_dh = dqsats(tsurf(i), zx_qs)/RLVTT &
                  / zx_pkh(i)
          ELSE
             zx_qs = qsatl(tsurf(i)) / ps(i)
             zx_dq_s_dh = dqsatl(tsurf(i), zx_qs)/RLVTT &
                  / zx_pkh(i)
          ENDIF
       ENDIF
       zx_dq_s_dt(i) = RCPD * zx_pkh(i) * zx_dq_s_dh
       zx_qsat(i) = zx_qs
       zx_coef(i) = coef1lay(i) &
            * (1.0+SQRT(u1lay(i)**2+v1lay(i)**2)) &
            * p1lay(i)/(RD*t1lay(i))

    ENDDO

    ! === Calcul de la temperature de surface ===

    ! zx_sl = chaleur latente d'evaporation ou de sublimation

    do i = 1, knon
       zx_sl(i) = RLVTT
       if (tsurf(i) .LT. RTT) zx_sl(i) = RLSTT
       zx_k1(i) = zx_coef(i)
    enddo

    do i = 1, knon
       ! Q
       zx_oq(i) = 1. - (beta(i) * zx_k1(i) * peqBcoef(i) * dtime)
       zx_mq(i) = beta(i) * zx_k1(i) * &
            (peqAcoef(i) - zx_qsat(i) &
            + zx_dq_s_dt(i) * tsurf(i)) &
            / zx_oq(i)
       zx_nq(i) = beta(i) * zx_k1(i) * (-1. * zx_dq_s_dt(i)) &
            / zx_oq(i)

       ! H
       zx_oh(i) = 1. - (zx_k1(i) * petBcoef(i) * dtime)
       zx_mh(i) = zx_k1(i) * petAcoef(i) / zx_oh(i)
       zx_nh(i) = - (zx_k1(i) * RCPD * zx_pkh(i))/ zx_oh(i)

       ! Tsurface
       tsurf_new(i) = (tsurf(i) + cal(i)/(RCPD * zx_pkh(i)) * dtime * &
            (radsol(i) + zx_mh(i) + zx_sl(i) * zx_mq(i)) &
            + dif_grnd(i) * t_grnd * dtime)/ &
            ( 1. - dtime * cal(i)/(RCPD * zx_pkh(i)) * ( &
            zx_nh(i) + zx_sl(i) * zx_nq(i)) &
            + dtime * dif_grnd(i))


       ! Y'a-t-il fonte de neige?

       ! fonte_neige = (nisurf /= is_oce) .AND. &
       ! & (snow(i) > epsfra .OR. nisurf == is_sic .OR. nisurf == is_lic) &
       ! & .AND. (tsurf_new(i) >= RTT)
       ! if (fonte_neige) tsurf_new(i) = RTT
       d_ts(i) = tsurf_new(i) - tsurf(i)
       ! zx_h_ts(i) = tsurf_new(i) * RCPD * zx_pkh(i)
       ! zx_q_0(i) = zx_qsat(i) + zx_dq_s_dt(i) * d_ts(i)
       !== flux_q est le flux de vapeur d'eau: kg/(m**2 s) positive vers bas
       !== flux_t est le flux de cpt (energie sensible): j/(m**2 s)
       evap(i) = - zx_mq(i) - zx_nq(i) * tsurf_new(i)
       fluxlat(i) = - evap(i) * zx_sl(i)
       fluxsens(i) = zx_mh(i) + zx_nh(i) * tsurf_new(i)
       ! Derives des flux dF/dTs (W m-2 K-1):
       dflux_s(i) = zx_nh(i)
       dflux_l(i) = (zx_sl(i) * zx_nq(i))
       ! Nouvelle valeure de l'humidite au dessus du sol
       qsat_new=zx_qsat(i) + zx_dq_s_dt(i) * d_ts(i)
       q1_new = peqAcoef(i) - peqBcoef(i)*evap(i)*dtime
       qsurf(i)=q1_new*(1.-beta(i)) + beta(i)*qsat_new
    ENDDO

  END SUBROUTINE calcul_fluxs

end module calcul_fluxs_m
