module time_phylmdz

  implicit none

  INTEGER, protected:: itap = 0 ! number of calls to "physiq"
  integer, protected:: itau_w ! pas de temps d'\'ecriture

contains

  subroutine increment_itap

    ! Incrémenter le compteur de la physique

    ! Libraries:
    use xios, only: xios_update_calendar

    use phyetat0_m, only: itau_phy

    !--------------------------------------------------

    itap = itap + 1
    itau_w = itau_phy + itap
    CALL xios_update_calendar(itap)

  end subroutine increment_itap

end module time_phylmdz
