subroutine conf_interface(tau_calv)

  ! From phylmd/conf_phys.F90,v 1.7 2005/07/05 07:21:23

  use IOIPSL
  implicit none

  ! Configuration de l'interace atm/surf
  !
  ! tau_calv:    temps de relaxation pour la fonte des glaciers

  REAL          :: tau_calv

  ! Local
  integer              :: numout = 6
  !
  !Config Key  = tau_calv
  !Config Desc = temps de relaxation pour fonte des glaciers en jours
  !Config Def  = 1 an 
  !Config Help = 
  !
  tau_calv = 360.*10.
  call getin('tau_calv',tau_calv)

  write(numout,*)' *********'
  WRITE(numout,*)' Configuration de l''interface atm/surfaces  : '
  WRITE(numout,*)' tau_calv = ',tau_calv
  return

end subroutine conf_interface
