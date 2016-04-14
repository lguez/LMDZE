module cv_thermo_m

  ! From LMDZ4/libf/phylmd/cvthermo.h, version 1.1.1.1 2004/05/19 12:53:09
  ! Thermodynamical constants for convectL

  implicit none

  real cpd, cpv, cl, rrv, rrd, lv0, g, rowl, t0
  real clmcpv, clmcpd, cpdmcp, cpvmcpd, cpvmcl
  real eps, epsi, epsim1
  real ginv, hrd
  real grav

contains

  SUBROUTINE cv_thermo

    ! Set thermodynamical constants for convectL

    use SUPHEC_M, only: rcpd, rcpv, rcw, rd, rg, rlvtt, rv

    !-----------------------------------------------------

    cpd = RCPD
    cpv = RCPV
    cl  = RCW
    rrv = RV
    rrd = RD
    lv0 = RLVTT
    g   = RG
    t0  = 273.15
    grav= g

    rowl=1000. ! (\`A quelle variable de SUPHEC_M cela correspond-il ?)
    clmcpv=cl-cpv
    clmcpd=cl-cpd
    cpdmcp=cpd-cpv
    cpvmcpd=cpv-cpd
    cpvmcl=cl-cpv
    eps=rrd/rrv
    epsi=1.0/eps
    epsim1=epsi-1.0
    ginv=1.0/grav
    hrd=0.5*rrd

  end SUBROUTINE cv_thermo

end module cv_thermo_m
