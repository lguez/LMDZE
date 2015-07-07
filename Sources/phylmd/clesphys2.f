module clesphys2

  ! From version 1.3 2005/06/06 13:16:33

  implicit none

  LOGICAL:: cycle_diurne = .TRUE.
  ! Cette option permet d'\'eteindre le cycle diurne. Peut \^etre utile
  ! pour acc\'el\'erer le code.

  LOGICAL:: soil_model = .TRUE.
  ! Choix du modele de sol (Thermique ?)

  LOGICAL:: new_oliq = .TRUE.
  ! Permet de mettre en route la nouvelle parametrisation de l'eau liquide

  ! Pour l'orographie:
  LOGICAL:: ok_orodr = .TRUE.
  LOGICAL:: ok_orolf = .TRUE.

  LOGICAL:: ok_limitvrai = .FALSE.
  ! On peut forcer le modele a lire le fichier SST de la bonne
  ! annee.

  INTEGER:: nbapp_rad = 12
  ! nombre d'appels des routines de rayonnements par jour

  INTEGER:: iflag_con = 2
  ! Convection scheme:
  ! 2 Tiedtke
  ! 3 Emanuel
  ! 4 Emanuel vect

contains

  subroutine read_clesphys2

    use unit_nml_m, only: unit_nml
    use nr_util, only: assert
    use conf_gcm_m, only: day_step, iphysiq

    namelist /clesphys2_nml/cycle_diurne, soil_model, new_oliq, ok_orodr, &
         ok_orolf, ok_limitvrai, nbapp_rad, iflag_con

    !------------------------------------

    print *, "Enter namelist 'clesphys2_nml'."
    read(unit=*, nml=clesphys2_nml)
    write(unit_nml, nml=clesphys2_nml)
    call assert(iflag_con >= 2 .and. iflag_con <= 4, "read_clesphys2 iflag_con")
    call assert(mod(day_step / iphysiq, nbapp_rad) == 0, &
         "read_clesphys2 nbapp_rad")
    call assert(nbapp_rad >= 4 .or. .not. cycle_diurne, &
         "read_clesphys2: minimum 4 calls to radiative transfer per day if " &
         // "cycle diurne")

  end subroutine read_clesphys2

end module clesphys2
