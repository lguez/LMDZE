module serre

  implicit none

  REAL:: clon = 0. ! longitude of the center of the zoom, in degrees
  real:: clat = 0. ! latitude of the center of the zoom, in degrees

  real:: grossismx = 1. ! facteur de grossissement du zoom, selon la longitude
  real:: grossismy = 1. ! facteur de grossissement du zoom, selon la latitude

  real:: dzoomx = 0.
  ! extension en longitude de la zone du zoom (fraction de la zone totale)

  real:: dzoomy = 0.
  ! extension en latitude de la zone du zoom (fraction de la zone totale)

  real:: taux = 3. ! raideur du zoom en X
  real:: tauy = 3. ! raideur du zoom en Y

contains

  subroutine read_serre

    use unit_nml_m, only: unit_nml
    use nr_util, only: assert

    namelist /serre_nml/ clon, clat, grossismx, grossismy, dzoomx, dzoomy, &
         taux, tauy

    !-------------------------------------------------

    print *, "Enter namelist 'serre_nml'."
    read(unit=*, nml=serre_nml)
    write(unit_nml, nml=serre_nml)

    call assert(grossismx >= 1. .and. grossismy >= 1., "read_serre grossism")

  end subroutine read_serre

end module serre
