
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/clift.F,v 1.1.1.1 2004/05/19
! 12:53:08 lmdzadmin Exp $

SUBROUTINE clift(p, t, rr, rs, plcl, dplcldt, dplcldq)
  ! ***************************************************************
  ! *                                                             *
  ! * CLIFT : COMPUTE LIFTING CONDENSATION LEVEL AND ITS          *
  ! *         DERIVATIVES RELATIVE TO T AND R                     *
  ! *   (WITHIN 0.2% OF FORMULA OF BOLTON, MON. WEA. REV.,1980)   *
  ! *                                                             *
  ! * written by   : GRANDPEIX Jean-Yves, 17/11/98, 12.39.01      *
  ! * modified by :                                               *
  ! ***************************************************************
  ! *
  ! *Arguments :
  ! *
  ! * Input :  P = pressure of level from wich lifting is performed
  ! *          T = temperature of level P
  ! *          RR = vapour mixing ratio at level P
  ! *          RS = vapour saturation mixing ratio at level P
  ! *
  ! * Output : PLCL = lifting condensation level
  ! *          DPLCLDT = derivative of PLCL relative to T
  ! *          DPLCLDQ = derivative of PLCL relative to R
  ! *
  ! cccccccccccccccccccccc
  ! constantes coherentes avec le modele du Centre Europeen
  ! RD = 1000.0 * 1.380658E-23 * 6.0221367E+23 / 28.9644
  ! RV = 1000.0 * 1.380658E-23 * 6.0221367E+23 / 18.0153
  ! CPD = 3.5 * RD
  ! CPV = 4.0 * RV
  ! CL = 4218.0
  ! CI=2090.0
  ! CPVMCL=CL-CPV
  ! CLMCI=CL-CI
  ! EPS=RD/RV
  ! ALV0=2.5008E+06
  ! ALF0=3.34E+05

  ! on utilise les constantes thermo du Centre Europeen: (sb)

  USE suphec_m

  cpd = rcpd
  cpv = rcpv
  cl = rcw
  cpvmcl = cl - cpv
  eps = rd/rv
  alv0 = rlvtt


  ! Bolton formula coefficients :
  a = 1669.0
  b = 122.0

  rh = rr/rs
  chi = t/(a-b*rh-t)
  plcl = p*(rh**chi)

  alv = alv0 - cpvmcl*(t-273.15)

  ! -- sb: correction:
  ! DPLCLDQ = PLCL*CHI*( 1./RR - B*CHI/T/RS*ALOG(RH) )
  dplcldq = plcl*chi*(1./rr+b*chi/t/rs*alog(rh))
  ! sb --

  dplcldt = plcl*chi*((a-b*rh*(1.+alv/rv/t))/t**2*chi*alog(rh)-alv/rv/t**2)


  RETURN
END SUBROUTINE clift
