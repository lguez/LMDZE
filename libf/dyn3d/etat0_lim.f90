PROGRAM etat0_lim

  ! This program sets the initial and boundary values.

  use comconst, only: initialize
  use conf_gcm_m, only: conf_gcm
  use etat0_mod, only: etat0
  use limit_mod, only: limit

  implicit none

  !-------------------------------------

  call initialize
  CALL conf_gcm
  CALL etat0
  CALL limit

END PROGRAM etat0_lim
