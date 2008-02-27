module iniprint

  ! gestion des impressions de sorties et de d�bogage

  implicit none

  INTEGER:: lunout = 6
  ! (unit� du fichier dans lequel se font les sorties, par d�faut la
  ! sortie standard)

  integer:: prt_level = 0 ! niveau d'impression souhait� (0 = minimum)

contains

  subroutine read_iniprint

    namelist /iniprint_nml/lunout, prt_level

    !-------------------------------------------------

    print *, "Call sequence information: read_iniprint"
    print *, "Enter namelist 'iniprint_nml'."
    read(unit=*, nml=iniprint_nml)
    write(unit=*, nml=iniprint_nml)

  end subroutine read_iniprint

end module iniprint
