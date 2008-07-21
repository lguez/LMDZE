module iniprint

  ! gestion des impressions de sorties et de débogage

  implicit none

  integer:: prt_level = 0 ! niveau d'impression souhaité (0 = minimum)

contains

  subroutine read_iniprint

    namelist /iniprint_nml/prt_level

    !-------------------------------------------------

    print *, "Call sequence information: read_iniprint"
    print *, "Enter namelist 'iniprint_nml'."
    read(unit=*, nml=iniprint_nml)
    write(unit=*, nml=iniprint_nml)

  end subroutine read_iniprint

end module iniprint
