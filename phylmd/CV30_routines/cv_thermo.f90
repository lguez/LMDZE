module cv_thermo

  ! From LMDZ4/libf/phylmd/cvthermo.h, version 1.1.1.1 2004/05/19 12:53:09
  ! Thermodynamical constants for cv_driver

  use SUPHEC_M, only: rd, rg, rcpv, rcw, rv

  implicit none

  real, parameter:: rowl = 1000.
  real, parameter:: clmcpv = rcw - rcpv
  real, parameter:: eps = rd / rv
  real, parameter:: ginv = 1. / rg

  private rd, rg, rcpv, rcw, rv

end module cv_thermo
