module time_phylmdz

  implicit none

  INTEGER:: itap = 0 ! number of calls to "physiq"
  integer itau_w ! pas de temps d'\'ecriture

contains

  subroutine increment_itap

    ! Incrémenter le compteur de la physique

    use phyetat0_m, only: itau_phy

    !--------------------------------------------------

    itap = itap + 1
    itau_w = itau_phy + itap

  end subroutine increment_itap

end module time_phylmdz
