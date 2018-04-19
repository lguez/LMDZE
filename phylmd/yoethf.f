module yoethf_m

  ! From phylmd/yoethf.inc, version 1.1.1.1 2004/05/19 12:53:09
  ! Derived constants specific to ECMWF thermodynamics

  use SUPHEC_M, only: ratm, rcpd, rcpv, rd, restt, rtt, rv

  implicit none

  ! Constants used for computation of saturation mixing ratio over
  ! liquid water (r.les) or ice (r.ies):
  real, parameter:: R3LES = 17.269, R3IES = 21.875, R4LES = 35.86, R4IES = 7.66
  real, parameter:: R5LES=R3LES*(RTT-R4LES)
  real, parameter:: R5IES=R3IES*(RTT-R4IES)

  real, parameter:: RVTMP2=RCPV/RCPD-1.
  real, parameter:: RHOH2O=RATM/100. ! density of liquid water
  real, parameter:: R2ES = RESTT * RD / RV

  private ratm, rcpd, rcpv, rd, restt, rtt, rv

end module yoethf_m
