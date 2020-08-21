module dimensions

  ! Model dimensions.

  implicit none

  INTEGER, protected:: iim = 16 ! number of longitudes
  INTEGER, protected:: jjm = 12 ! number of latitudes

  INTEGER, protected:: llm = 11
  ! number of vertical layers, should be >= 7, because of nl in cv30_param_m

  integer, parameter:: nqmx = 5 ! maximum number of tracers

contains

  subroutine set_dimensions

    use unit_nml_m, only: unit_nml

    namelist /dimensions_nml/iim, jjm, llm

    !------------------------------------------------------------

    print *, "Enter namelist 'dimensions_nml'."
    read (unit=*, nml=dimensions_nml)
    write(unit_nml, nml=dimensions_nml)

  end subroutine set_dimensions
  
end module dimensions
