module conf_guide_m

  IMPLICIT NONE

  LOGICAL:: ini_anal = .false. ! Initial = analyse
  LOGICAL:: guide_u = .false. ! guidage de u
  LOGICAL:: guide_v = .false. ! gvidage de v
  LOGICAL:: guide_t = .false. ! guidage de T
  LOGICAL:: guide_q = .false. ! guidage de q

  logical, save:: ok_guide ! guidage

  ! alpha détermine la part des injections de données à chaque étape
  ! alpha=0 signifie pas d'injection
  ! alpha=1 signifie injection totale
  REAL, save, allocatable, protected:: alpha_q(:, :), alpha_t(:, :), &
       alpha_u(:, :) ! (iim + 1, jjm + 1)
  REAL, save, allocatable, protected:: alpha_v(:, :) ! (iim + 1, jjm)

contains

  SUBROUTINE conf_guide

    ! From LMDZ4/libf/dyn3d/conf_guide.F, version 1.1.1.1 2004/05/19 12:53:07
    !  Parametres de controle du run:

    use comconst, only: daysec
    USE dimensions, ONLY: iim, jjm
    use conf_gcm_m, only: day_step, iperiod, dtvr
    use dynetat0_m, only: rlatu, rlatv
    use dynetat0_chosen_m, only: grossismx, grossismy
    use init_tau2alpha_m, only: init_tau2alpha
    use jumble, only: assert, pi
    use tau2alpha_m, only: tau2alpha
    use unit_nml_m, only: unit_nml
    use writefield_m, only: writefield

    ! Local:

    ! Constantes de rappel, en jours :
    REAL tau_min_u, tau_max_u, tau_min_v, tau_max_v, tau_min_t, tau_max_t
    real tau_min_q, tau_max_q

    REAL lat_min_guide_deg, lat_max_guide_deg ! in degrees

    ! Dans le cas où on n'a les analyses que sur une bande de latitudes :
    REAL lat_min_guide ! minimum latitude for nudging, in rad
    real lat_max_guide ! maximum latitude for nudging, in rad

    REAL factt ! pas de temps entre deux appels au guidage, en jours
    REAL dxdys(iim + 1, jjm + 1), dxdyu(iim + 1, jjm + 1), dxdyv(iim + 1, jjm)
    integer i, j

    namelist /conf_guide_nml/ ini_anal, guide_u, guide_v, guide_t, guide_q, &
         tau_min_u, tau_max_u, tau_min_v, tau_max_v, tau_min_t, tau_max_t, &
         tau_min_q, tau_max_q, lat_min_guide_deg, lat_max_guide_deg

    !-----------------------------------------------------------------------

    print *, "Call sequence information: conf_guide"

    allocate(alpha_q(iim + 1, jjm + 1))
    allocate(alpha_t(iim + 1, jjm + 1))
    allocate(alpha_u(iim + 1, jjm + 1), alpha_v(iim + 1, jjm))

    ! Default values:
    tau_min_u = 0.03
    tau_max_u = 10.
    tau_min_v = 0.03
    tau_max_v = 10.
    tau_min_t = 0.03
    tau_max_t = 10.
    tau_min_q = 0.03
    tau_max_q = 10.
    lat_min_guide_deg = -90.
    lat_max_guide_deg = 90.

    print *, "Enter namelist 'conf_guide_nml'."
    read(unit=*, nml=conf_guide_nml)
    write(unit_nml, nml=conf_guide_nml)

    ok_guide = any((/guide_u, guide_v, guide_t, guide_q/))

    if (ok_guide) then
       call assert(mod(day_step, 4 * iperiod) == 0, &
            "conf_guide ok_guide day_step iperiod")
       lat_min_guide = lat_min_guide_deg / 180. * pi
       lat_max_guide = lat_max_guide_deg / 180. * pi
       factt = dtvr * iperiod / daysec
       print *, "factt = ", factt

       IF (abs(grossismx - 1.) < 0.1 .OR. abs(grossismy - 1.) < 0.1) THEN
          ! grille regulière
          if (guide_u) alpha_u = 1. - exp(- factt / tau_max_u)
          if (guide_v) alpha_v = 1. - exp(- factt / tau_max_v)
          if (guide_t) alpha_t = 1. - exp(- factt / tau_max_t)
          if (guide_q) alpha_q = 1. - exp(- factt / tau_max_q)
       else
          call init_tau2alpha(dxdys)

          if (guide_u) then
             DO j = 1, jjm + 1
                DO i = 1, iim
                   dxdyu(i, j) = 0.5 * (dxdys(i, j) + dxdys(i + 1, j))
                END DO
                dxdyu(iim + 1, j) = dxdyu(1, j)
             END DO

             CALL tau2alpha(lat_min_guide, lat_max_guide, factt, dxdyu, rlatu, &
                  tau_min_u, tau_max_u, alpha_u)
             CALL writefield("alpha_u", alpha_u)
          end if

          if (guide_v) then
             DO j = 1, jjm
                DO i = 1, iim + 1
                   dxdyv(i, j) = 0.5 * (dxdys(i, j) + dxdys(i, j + 1))
                END DO
             END DO

             CALL tau2alpha(lat_min_guide, lat_max_guide, factt, dxdyv, rlatv, &
                  tau_min_v, tau_max_v, alpha_v)
             CALL writefield("alpha_v", alpha_v)
          end if

          if (guide_t) then
             CALL tau2alpha(lat_min_guide, lat_max_guide, factt, dxdys, rlatu, &
                  tau_min_t, tau_max_t, alpha_t)
             CALL writefield("alpha_t", alpha_t)
          end if

          if (guide_q)  then
             CALL tau2alpha(lat_min_guide, lat_max_guide, factt, dxdys, rlatu, &
                  tau_min_q, tau_max_q, alpha_q)
             CALL writefield("alpha_q", alpha_q)
          end if
       end IF
    end if

  end SUBROUTINE conf_guide

end module conf_guide_m
