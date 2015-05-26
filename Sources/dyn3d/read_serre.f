module read_serre_m

  implicit none

contains

  subroutine read_serre

    use dynetat0_m, only: clon, clat, grossismx, grossismy, dzoomx, dzoomy, &
         taux, tauy
    use unit_nml_m, only: unit_nml
    use nr_util, only: assert, pi

    REAL:: clon_deg = 0. ! longitude of the center of the zoom, in degrees
    real:: clat_deg = 0. ! latitude of the center of the zoom, in degrees

    namelist /serre_nml/ clon_deg, clat_deg, grossismx, grossismy, dzoomx, &
         dzoomy, taux, tauy

    !-------------------------------------------------

    ! Default values:
    clon_deg = 0.
    clat_deg = 0.
    grossismx = 1.
    grossismy = 1.
    dzoomx = 0.2
    dzoomy = 0.2
    taux = 3.
    tauy = 3. 

    print *, "Enter namelist 'serre_nml'."
    read(unit=*, nml=serre_nml)
    write(unit_nml, nml=serre_nml)

    call assert(grossismx >= 1. .and. grossismy >= 1., "read_serre grossism")
    call assert(dzoomx > 0., dzoomx < 1., dzoomy < 1., &
         "read_serre dzoomx dzoomy")
    clon = clon_deg / 180. * pi
    clat = clat_deg / 180. * pi

  end subroutine read_serre

end module read_serre_m
