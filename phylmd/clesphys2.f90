module clesphys2

  ! From version 1.3 2005/06/06 13:16:33

  implicit none

  LOGICAL, protected:: new_oliq = .TRUE.
  ! Permet de mettre en route la nouvelle parametrisation de l'eau liquide

  ! Pour l'orographie:
  LOGICAL, protected:: ok_orodr = .TRUE.
  LOGICAL, protected:: ok_orolf = .TRUE.

  LOGICAL, protected:: ok_limitvrai = .FALSE.
  ! On peut forcer le modele a lire le fichier SST de la bonne
  ! annee.

  INTEGER, protected:: nbapp_rad = 12
  ! nombre d'appels des routines de rayonnements par jour

  logical, protected:: conv_emanuel = .true.
  ! convection scheme of Emanuel, else Tiedtke

contains

  subroutine read_clesphys2

    use unit_nml_m, only: unit_nml
    use nr_util, only: assert
    use conf_gcm_m, only: day_step, iphysiq

    namelist /clesphys2_nml/new_oliq, ok_orodr, ok_orolf, ok_limitvrai, &
         nbapp_rad, conv_emanuel

    !------------------------------------

    print *, "Enter namelist 'clesphys2_nml'."
    read(unit=*, nml=clesphys2_nml)
    write(unit_nml, nml=clesphys2_nml)
    call assert(mod(day_step / iphysiq, nbapp_rad) == 0, &
         "read_clesphys2 nbapp_rad")
    call assert(nbapp_rad >= 4, &
         "read_clesphys2: minimum 4 calls to radiative transfer per day")

  end subroutine read_clesphys2

end module clesphys2
