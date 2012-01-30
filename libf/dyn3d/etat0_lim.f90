PROGRAM etat0_lim

  ! This program sets the initial and boundary values.

  use conf_gcm_m, only: conf_gcm
  use etat0_mod, only: etat0
  use jumble, only: new_unit
  use limit_mod, only: limit
  use unit_nml_m, only: unit_nml

  implicit none

  !-------------------------------------

  call new_unit(unit_nml)
  open(unit_nml, file="used_namelists", status="replace", action="write")
  CALL conf_gcm
  CALL etat0
  CALL limit
  close(unit_nml)

END PROGRAM etat0_lim
