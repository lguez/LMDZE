SUBROUTINE cv_thermo

  ! Set thermodynamical constants for convectL

  use clesphys2, only: iflag_con
  use SUPHEC_M
  use cvthermo

  implicit none

  !-----------------------------------------------------

  ! original set from convect:
  if (iflag_con.eq.4) then
     cpd=1005.7
     cpv=1870.0
     cl=4190.0
     rrv=461.5
     rrd=287.04
     lv0=2.501E6
     g=9.8
     t0=273.15
     grav=g
  endif

  ! constants consistent with LMDZ:
  if (iflag_con.eq.3) then
     cpd = RCPD
     cpv = RCPV
     cl  = RCW
     rrv = RV
     rrd = RD
     lv0 = RLVTT
     g   = RG     ! not used in convect3
     t0  = 273.15 ! convect3 (RTT=273.16)
     grav= g    ! implicitely or explicitely used in convect3
  endif

  rowl=1000.0 !(a quelle variable de SUPHEC_M cela correspond-il?)
  clmcpv=cl-cpv
  clmcpd=cl-cpd
  cpdmcp=cpd-cpv
  cpvmcpd=cpv-cpd
  cpvmcl=cl-cpv ! for convect3
  eps=rrd/rrv
  epsi=1.0/eps
  epsim1=epsi-1.0
  ginv=1.0/grav
  hrd=0.5*rrd

end SUBROUTINE cv_thermo
