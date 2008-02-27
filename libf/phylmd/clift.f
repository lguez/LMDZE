!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/clift.F,v 1.1.1.1 2004/05/19 12:53:08 lmdzadmin Exp $
!
        SUBROUTINE CLIFT (P,T,RR,RS,PLCL,DPLCLDT,DPLCLDQ)
C***************************************************************
C*                                                             *
C* CLIFT : COMPUTE LIFTING CONDENSATION LEVEL AND ITS          *
C*         DERIVATIVES RELATIVE TO T AND R                     *
C*   (WITHIN 0.2% OF FORMULA OF BOLTON, MON. WEA. REV.,1980)   *
C*                                                             *
C* written by   : GRANDPEIX Jean-Yves, 17/11/98, 12.39.01      *
C* modified by :                                               *
C***************************************************************
C*
C*Arguments :
C*
C* Input :  P = pressure of level from wich lifting is performed
C*          T = temperature of level P
C*          RR = vapour mixing ratio at level P
C*          RS = vapour saturation mixing ratio at level P
C*
C* Output : PLCL = lifting condensation level
C*          DPLCLDT = derivative of PLCL relative to T
C*          DPLCLDQ = derivative of PLCL relative to R
C*
ccccccccccccccccccccccc
c constantes coherentes avec le modele du Centre Europeen
c      RD = 1000.0 * 1.380658E-23 * 6.0221367E+23 / 28.9644
c      RV = 1000.0 * 1.380658E-23 * 6.0221367E+23 / 18.0153
c      CPD = 3.5 * RD
c      CPV = 4.0 * RV
c      CL = 4218.0
c      CI=2090.0
c      CPVMCL=CL-CPV
c      CLMCI=CL-CI
c      EPS=RD/RV
c      ALV0=2.5008E+06
c      ALF0=3.34E+05
c
c on utilise les constantes thermo du Centre Europeen: (sb)
c
      use YOMCST
c
       CPD = RCPD
       CPV = RCPV
       CL = RCW
       CPVMCL = CL-CPV
       EPS = RD/RV
       ALV0 = RLVTT
c
c
c      Bolton formula coefficients :
      A = 1669.0
      B = 122.0
c
      RH=RR/RS
      CHI=T/(A-B*RH-T)
      PLCL=P*(RH**CHI)
c
      ALV = ALV0 - CPVMCL*(T-273.15)
c
c -- sb: correction:
c       DPLCLDQ = PLCL*CHI*( 1./RR - B*CHI/T/RS*ALOG(RH) )
      DPLCLDQ = PLCL*CHI*( 1./RR + B*CHI/T/RS*ALOG(RH) )
c sb --
c
      DPLCLDT = PLCL*CHI*((A-B*RH*(1.+ALV/RV/T))/T**2*CHI*ALOG(RH)
     $                    - ALV/RV/T**2 )
c
c
      RETURN
      END
