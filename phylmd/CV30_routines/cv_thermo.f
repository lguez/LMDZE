module cv_thermo

  ! From LMDZ4/libf/phylmd/cvthermo.h, version 1.1.1.1 2004/05/19 12:53:09
  ! Thermodynamical constants for cv_driver

  use SUPHEC_M, only: rd, rg, rcpd, rcpv, rcw, rv

  implicit none

  real, parameter:: rowl = 1000., t0 = 273.15
  real, parameter:: clmcpv = rcw - rcpv, clmcpd = rcw - rcpd
  real, parameter:: cpdmcp = rcpd - rcpv
  real, parameter:: cpvmcpd = rcpv - rcpd, cpvmcl = rcw - rcpv
  real, parameter:: eps = rd / rv, epsi = 1. / eps, epsim1 = epsi - 1.
  real, parameter:: ginv = 1. / rg, hrd = 0.5 * rd

  private rd, rg, rcpd, rcpv, rcw, rv

end module cv_thermo
