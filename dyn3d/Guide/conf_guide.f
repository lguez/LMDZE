module conf_guide_m

  IMPLICIT NONE

  !   Constantes de rappel. Unite : fraction de jour
  REAL:: tau_min_u = 0.02
  REAL:: tau_max_u = 10.
  REAL:: tau_min_v = 0.02
  REAL:: tau_max_v = 10.
  REAL:: tau_min_t = 0.02
  REAL:: tau_max_t = 10.
  REAL:: tau_min_q = 0.02
  REAL:: tau_max_q = 10.
  REAL:: tau_min_p = 0.02
  REAL:: tau_max_p = 10.

  LOGICAL:: ncep = .false. ! Coordonnee vert NCEP ou ECMWF
  LOGICAL:: ini_anal = .false. ! Initial = analyse
  LOGICAL:: guide_u = .true. ! guidage de u
  LOGICAL:: guide_v = .true. ! gvidage de v
  LOGICAL:: guide_t = .false. ! guidage de T
  LOGICAL:: guide_q = .false. ! guidage de q

  logical:: online = .true. ! controle du guide
  ! hors-ligne: x=x_rea

  ! Latitude min et max pour le rappel dans le cas ou on 'a les
  ! analyses que sur une bande de latitudes.
  REAL:: lat_min_guide = -90. ! Latitude minimum pour le guidage
  real:: lat_max_guide = 90. ! Latitude maximum pour le guidage

contains

  SUBROUTINE conf_guide

    ! From LMDZ4/libf/dyn3d/conf_guide.F, version 1.1.1.1 2004/05/19 12:53:07
    !  Parametres de controle du run:

    use unit_nml_m, only: unit_nml

    namelist /conf_guide_nml/ ncep, ini_anal, guide_u, guide_v, guide_t, &
         online, tau_min_u, tau_max_u, tau_min_v, tau_max_v, tau_min_t, &
         tau_max_t, tau_min_q, tau_max_q, tau_min_p, tau_max_p, &
         lat_min_guide, lat_max_guide

    !-----------------------------------------------------------------------

    print *, "Call sequence information: conf_guide"

    print *, "Enter namelist 'conf_guide_nml'."
    read(unit=*, nml=conf_guide_nml)
    write(unit_nml, nml=conf_guide_nml)

  end SUBROUTINE conf_guide

end module conf_guide_m
