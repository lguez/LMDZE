SUBROUTINE cv_thermo

  ! Set thermodynamical constants for convectL

  use SUPHEC_M
  use cvthermo

  implicit none

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

  rowl=1000.0 !(a quelle variable de SUPHEC_M cela correspond-il?)
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
