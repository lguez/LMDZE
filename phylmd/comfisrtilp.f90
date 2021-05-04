module comfisrtilp

  ! From phylmd/fisrtilp.h,v 1.1.1.1 2004/05/19 12:53:09

  implicit none

  REAL, protected:: cld_lc_lsc  = 2.6e-4
  real, protected:: cld_lc_con = 2.6e-4
  real, protected:: cld_tau_lsc = 3600.
  real, protected:: cld_tau_con = 3600.
  real, protected:: ffallv_lsc = 1.
  real, protected:: ffallv_con = 1.
  real, protected:: coef_eva = 2.e-5
  logical, protected:: reevap_ice = .false.
  integer, protected:: iflag_pdf = 0

contains

  subroutine read_comfisrtilp

    use unit_nml_m, only: unit_nml

    namelist /fisrtilp_nml/ cld_lc_lsc, cld_lc_con, cld_tau_lsc, cld_tau_con, &
         ffallv_lsc, ffallv_con, coef_eva, reevap_ice, iflag_pdf

    !----------------------------------------------------------------

    print *, "Enter namelist 'fisrtilp_nml'."
    read(unit=*, nml = fisrtilp_nml)
    write(unit_nml, nml = fisrtilp_nml)

  end subroutine read_comfisrtilp

end module comfisrtilp
