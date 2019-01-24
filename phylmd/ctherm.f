module ctherm_m

  implicit none

  logical, protected:: iflag_thermals = .false.
  integer, protected:: nsplit_thermals = 1

contains

  subroutine ctherm

    use unit_nml_m, only: unit_nml

    namelist /ctherm_nml/ iflag_thermals, nsplit_thermals

    !--------------------------------------------------------------------

    print *, "Enter namelist 'ctherm_nml'."
    read(unit=*, nml=ctherm_nml)
    write(unit_nml, nml=ctherm_nml)

  end subroutine ctherm

end module ctherm_m
