module serre

  implicit none

  REAL:: clon = 0. ! longitude of the center of the zoom, in degrees
  real:: clat = 0. ! latitude of the center of the zoom, in degrees

  real:: grossismx = 1., grossismy = 1.
  ! facteurs de grossissement du zoom, selon la longitude et la latitude
  ! = 2 si 2 fois, = 3 si 3 fois, etc.

  real:: dzoomx = 0.2, dzoomy = 0.2
  ! extensions en longitude et latitude de la zone du zoom (fractions
  ! de la zone totale)

  real:: taux = 3., tauy = 3. 
  ! raideur de la transition de l'intérieur à l'extérieur du zoom
  
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
    call assert(dzoomx > 0., dzoomx < 1., dzoomy < 1., &
         "read_serre dzoomx dzoomy")

  end subroutine read_serre

end module serre
