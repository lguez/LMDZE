module conf_phys_m

  implicit none

  integer, protected:: iflag_pbl = 1 ! for the planetary boundary layer
  ! Allowed values: 0, 1, 6, 8, 9
  ! 6 : Mellor and Yamada 2.0
  ! 8 : Mellor and Yamada 2.5

  REAL, protected:: rad_chau1 = 13., rad_chau2 = 9.
  ! effective radii of liquid cloud droplets, in micrometers
  
  real:: epmax = 0.993 ! \'efficacit\'e de pr\'ecipitation
  integer:: iflag_clw = 0
  integer, protected, save:: kdlon

contains

  subroutine conf_phys

    ! From phylmd/conf_phys.F90, version 1.7 2005/07/05 07:21:23

    ! Configuration de la "physique" de LMDZ.

    USE clesphys, ONLY: read_clesphys
    use clesphys2, only: read_clesphys2
    use comfisrtilp, only: read_comfisrtilp
    use dimphy, only: klon
    use nr_util, only: assert
    use unit_nml_m, only: unit_nml
    USE yomcst, ONLY: read_YOMCST

    namelist /conf_phys_nml/ epmax, iflag_clw, iflag_pbl

    namelist /nuagecom/ rad_chau1, rad_chau2

    !-----------------------------------------------------------

    print *, "Call sequence information: conf_phys"
    kdlon = klon
    call read_clesphys2
    call read_YOMCST


    print *, "Enter namelist 'conf_phys_nml'."
    read(unit=*, nml=conf_phys_nml)
    write(unit_nml, nml=conf_phys_nml)

    call assert(any(iflag_pbl == [0, 1, 6, 8, 9]), &
         "conf_phys: bad value for iflag_pbl")

    call read_comfisrtilp
    call read_clesphys

    print *, "Enter namelist 'nuagecom'."
    read(unit=*, nml=nuagecom)
    write(unit_nml, nml=nuagecom)

  end subroutine conf_phys

end module conf_phys_m
