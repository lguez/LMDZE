module conf_phys_m

  implicit none

  integer, protected:: iflag_pbl = 1 ! for the planetary boundary layer
  ! 6 : Mellor and Yamada 2.0
  ! 8 : Mellor and Yamada 2.5

  REAL, protected:: rad_chau1 = 13., rad_chau2 = 9.
  real:: epmax = 0.993 ! \'efficacit\'e de pr\'ecipitation
  integer:: iflag_clw = 0

contains

  subroutine conf_phys

    ! From phylmd/conf_phys.F90, version 1.7 2005/07/05 07:21:23

    ! Configuration de la "physique" de LMDZ.

    USE clesphys, ONLY: read_clesphys
    use clesphys2, only: read_clesphys2
    USE comfisrtilp, ONLY: cld_lc_con, cld_lc_lsc, cld_tau_con, &
         cld_tau_lsc, coef_eva, ffallv_con, ffallv_lsc, iflag_pdf, reevap_ice
    use nr_util, only: assert
    use unit_nml_m, only: unit_nml
    USE yomcst, ONLY: read_YOMCST

    namelist /conf_phys_nml/ epmax, iflag_clw, cld_lc_lsc, cld_lc_con, &
         cld_tau_lsc, cld_tau_con, ffallv_lsc, ffallv_con, coef_eva, &
         reevap_ice, iflag_pdf, iflag_pbl

    namelist /nuagecom/ rad_chau1, rad_chau2

    !-----------------------------------------------------------

    print *, "Call sequence information: conf_phys"
    call read_clesphys2
    call read_YOMCST

    cld_lc_lsc = 2.6e-4
    cld_lc_con = 2.6e-4
    cld_tau_lsc = 3600.
    cld_tau_con = 3600.
    ffallv_lsc = 1.
    ffallv_con = 1.
    coef_eva = 2.e-5
    reevap_ice = .false.
    iflag_pdf = 0

    print *, "Enter namelist 'conf_phys_nml'."
    read(unit=*, nml=conf_phys_nml)
    write(unit_nml, nml=conf_phys_nml)

    call assert(any(iflag_pbl == [0, 1, 6, 8, 9]), &
         "conf_phys: bad value for iflag_pbl")
    call read_clesphys

    print *, "Enter namelist 'nuagecom'."
    read(unit=*, nml=nuagecom)
    write(unit_nml, nml=nuagecom)

  end subroutine conf_phys

end module conf_phys_m
