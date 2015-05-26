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
  REAL:: tau_min_p = 0.03
  REAL:: tau_max_p = 10.

  LOGICAL:: ncep = .false. ! Coordonnee vert NCEP ou ECMWF
  LOGICAL:: ini_anal = .false. ! Initial = analyse
  LOGICAL:: guide_u = .false. ! guidage de u
  LOGICAL:: guide_v = .false. ! gvidage de v
  LOGICAL:: guide_t = .false. ! guidage de T
  LOGICAL:: guide_q = .false. ! guidage de q

  logical:: online = .true. ! controle du guide
  ! hors-ligne: x=x_rea

  ! Latitude min et max pour le rappel dans le cas ou on 'a les
  ! analyses que sur une bande de latitudes.
  REAL:: lat_min_guide = -90. ! Latitude minimum pour le guidage
  real:: lat_max_guide = 90. ! Latitude maximum pour le guidage

  logical, save:: ok_guide ! guidage
  REAL, save:: factt ! pas de temps entre deux appels au guidage, en jours

contains

  SUBROUTINE conf_guide

    ! From LMDZ4/libf/dyn3d/conf_guide.F, version 1.1.1.1 2004/05/19 12:53:07
    !  Parametres de controle du run:

    use abort_gcm_m, only: abort_gcm
    use comconst, only: daysec, dtvr
    use conf_gcm_m, only: day_step, iperiod
    use dynetat0_m, only: grossismx, grossismy
    use nr_util, only: assert
    use unit_nml_m, only: unit_nml

    namelist /conf_guide_nml/ ncep, ini_anal, guide_u, guide_v, guide_t, &
         guide_q, online, tau_min_u, tau_max_u, tau_min_v, tau_max_v, &
         tau_min_t, tau_max_t, tau_min_q, tau_max_q, tau_min_p, tau_max_p, &
         lat_min_guide, lat_max_guide

    !-----------------------------------------------------------------------

    print *, "Call sequence information: conf_guide"

    print *, "Enter namelist 'conf_guide_nml'."
    read(unit=*, nml=conf_guide_nml)
    write(unit_nml, nml=conf_guide_nml)

    ok_guide = any((/guide_u, guide_v, guide_t, guide_q/))
    if (ok_guide .and. mod(day_step, 4 * iperiod) /= 0) call &
         abort_gcm(modname = "conf_guide", &
         message = 'ok_guide day_step iperiod', ierr = 1)

    if (ok_guide .and. online) then
       factt = dtvr * iperiod / daysec
       print *, "factt = ", factt
       if (abs(grossismx - 1.) >= 0.1 .and. abs(grossismy - 1.) >= 0.1) then
          if (guide_u) call assert(factt / tau_min_u < 1, &
               "conf_guide tau_min_u")
          if (guide_v) call assert(factt / tau_min_v < 1, &
               "conf_guide tau_min_v")
          if (guide_t) call assert(factt / tau_min_t < 1, &
               "conf_guide tau_min_t")
          if (guide_q) call assert(factt / tau_min_q < 1, &
               "conf_guide tau_min_q")
       end if
    end if

  end SUBROUTINE conf_guide

end module conf_guide_m
