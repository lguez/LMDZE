module cv_thermo_m

  ! From LMDZ4/libf/phylmd/cvthermo.h, version 1.1.1.1 2004/05/19 12:53:09
  ! Thermodynamical constants for cv_driver

  use SUPHEC_M, only: rd, rg, rcpd, rcpv, rcw

  implicit none

  real cpd, cpv, cl, rrv, rrd, rowl, t0
  real, parameter:: clmcpv = rcw - rcpv, clmcpd = rcw - rcpd
  real, parameter:: cpdmcp = rcpd - rcpv
  real, parameter:: cpvmcpd = rcpv - rcpd, cpvmcl = rcw - rcpv
  real eps, epsi, epsim1
  real, parameter:: ginv = 1. / rg, hrd = 0.5 * rd

  private rd, rg, rcpd, rcpv, rcw

contains

  SUBROUTINE cv_thermo

    ! Set thermodynamical constants for cv_driver

    use SUPHEC_M, only: rlvtt, rv

    !-----------------------------------------------------

    cpd = RCPD
    cpv = RCPV
    cl  = RCW
    rrv = RV
    rrd = RD
    t0  = 273.15

    rowl = 1000. ! (\`A quelle variable de SUPHEC_M cela correspond-il ?)
    eps = rd/rrv
    epsi = 1.0/eps
    epsim1 = epsi - 1.0

  end SUBROUTINE cv_thermo

end module cv_thermo_m
