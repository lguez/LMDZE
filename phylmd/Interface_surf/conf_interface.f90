MODULE conf_interface_m

  ! From phylmd/interface_surf.F90, version 1.8 2005/05/25 13:10:09
  ! L. Fairhead, LMD, february 2000

  IMPLICIT none

  REAL, save, protected:: tau_calv 
  ! temps de relaxation pour la fonte des glaciers, en jours

contains

  subroutine conf_interface

    ! From phylmd/conf_phys.F90, version 1.7 2005/07/05 07:21:23

    ! Configuration de l'interace atm/surf

    use unit_nml_m, only: unit_nml

    namelist /conf_interface_nml/ tau_calv

    !------------------------------------------------------

    tau_calv = 360.*10.

    print *, "Enter namelist 'conf_interface_nml'."
    read(unit=*, nml=conf_interface_nml)
    write(unit_nml, nml=conf_interface_nml)

  end subroutine conf_interface

END MODULE conf_interface_m
