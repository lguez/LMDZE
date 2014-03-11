module yoethf_m

  ! From phylmd/yoethf.inc, version 1.1.1.1 2004/05/19 12:53:09
  ! Derived constants specific to ECMWF thermodynamics

  implicit none

  ! Constants used for computation of saturation mixing ratio over
  ! liquid water (r.les) or ice (r.ies):
  REAL R2ES, R5LES, R5IES
  real, parameter:: R3LES = 17.269, R3IES = 21.875, R4LES = 35.86, R4IES = 7.66

  REAL RVTMP2
  real RHOH2O ! density of liquid water

  save

contains

  subroutine yoethf

    ! Calcul des constantes pour les fonctions thermodynamiques

    use SUPHEC_M, only: ratm, rcpd, rcpv, rd, restt, rtt, rv

    !-----------------------------------------

    RVTMP2=RCPV/RCPD-1.
    RHOH2O=RATM/100.
    R2ES = RESTT * RD / RV
    R5LES=R3LES*(RTT-R4LES)
    R5IES=R3IES*(RTT-R4IES)

  end subroutine yoethf

end module yoethf_m
