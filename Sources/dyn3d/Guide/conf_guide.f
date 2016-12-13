module conf_guide_m

  IMPLICIT NONE

  ! Constantes de rappel, en jours :
  REAL:: tau_min_u = 0.03
  REAL:: tau_max_u = 10.
  REAL:: tau_min_v = 0.03
  REAL:: tau_max_v = 10.
  REAL:: tau_min_t = 0.03
  REAL:: tau_max_t = 10.
  REAL:: tau_min_q = 0.03
  REAL:: tau_max_q = 10.

  LOGICAL:: ini_anal = .false. ! Initial = analyse
  LOGICAL:: guide_u = .false. ! guidage de u
  LOGICAL:: guide_v = .false. ! gvidage de v
  LOGICAL:: guide_t = .false. ! guidage de T
  LOGICAL:: guide_q = .false. ! guidage de q

  logical:: online = .true. ! contr\^ole du guidage
  ! hors-ligne: x = x_rea

  ! Dans le cas où on n'a les analyses que sur une bande de latitudes :
  REAL, save:: lat_min_guide ! minimum latitude for nudging, in rad
  real, save:: lat_max_guide ! maximum latitude for nudging, in rad

  logical, save:: ok_guide ! guidage
  REAL, save:: factt ! pas de temps entre deux appels au guidage, en jours

contains

  SUBROUTINE conf_guide

    ! From LMDZ4/libf/dyn3d/conf_guide.F, version 1.1.1.1 2004/05/19 12:53:07
    !  Parametres de controle du run:

    use abort_gcm_m, only: abort_gcm
    use comconst, only: daysec, dtvr
    use conf_gcm_m, only: day_step, iperiod
    use nr_util, only: assert, pi
    use unit_nml_m, only: unit_nml

    ! Local:

    REAL:: lat_min_guide_deg = -90. ! in degrees
    real:: lat_max_guide_deg = 90. ! in degrees

    namelist /conf_guide_nml/ ini_anal, guide_u, guide_v, guide_t, guide_q, &
         online, tau_min_u, tau_max_u, tau_min_v, tau_max_v, tau_min_t, &
         tau_max_t, tau_min_q, tau_max_q, lat_min_guide_deg, lat_max_guide_deg

    !-----------------------------------------------------------------------

    print *, "Call sequence information: conf_guide"

    print *, "Enter namelist 'conf_guide_nml'."
    read(unit=*, nml=conf_guide_nml)
    write(unit_nml, nml=conf_guide_nml)

    lat_min_guide = lat_min_guide_deg / 180. * pi
    lat_max_guide = lat_max_guide_deg / 180. * pi

    ok_guide = any((/guide_u, guide_v, guide_t, guide_q/))
    if (ok_guide .and. mod(day_step, 4 * iperiod) /= 0) call &
         abort_gcm("conf_guide", 'ok_guide day_step iperiod')

    if (ok_guide .and. online) then
       factt = dtvr * iperiod / daysec
       print *, "factt = ", factt
    end if

  end SUBROUTINE conf_guide

end module conf_guide_m
