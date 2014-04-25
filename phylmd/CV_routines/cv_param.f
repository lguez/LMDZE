!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/cv_routines.F,v 1.1.1.1 2004/05/19 12:53:08 lmdzadmin Exp $
!
      SUBROUTINE cv_param(nd)
            use cvparam
      implicit none

!------------------------------------------------------------
! Set parameters for convectL
! (includes microphysical parameters and parameters that
!  control the rate of approach to quasi-equilibrium)
!------------------------------------------------------------

!   *** ELCRIT IS THE AUTOCONVERSION THERSHOLD WATER CONTENT (gm/gm) ***
!   ***  TLCRIT IS CRITICAL TEMPERATURE BELOW WHICH THE AUTO-        ***
!   ***       CONVERSION THRESHOLD IS ASSUMED TO BE ZERO             ***
!   ***     (THE AUTOCONVERSION THRESHOLD VARIES LINEARLY            ***
!   ***               BETWEEN 0 C AND TLCRIT)                        ***
!   ***   ENTP IS THE COEFFICIENT OF MIXING IN THE ENTRAINMENT       ***
!   ***                       FORMULATION                            ***
!   ***  SIGD IS THE FRACTIONAL AREA COVERED BY UNSATURATED DNDRAFT  ***
!   ***  SIGS IS THE FRACTION OF PRECIPITATION FALLING OUTSIDE       ***
!   ***                        OF CLOUD                              ***
!   ***        OMTRAIN IS THE ASSUMED FALL SPEED (P/s) OF RAIN       ***
!   ***     OMTSNOW IS THE ASSUMED FALL SPEED (P/s) OF SNOW          ***
!   ***  COEFFR IS A COEFFICIENT GOVERNING THE RATE OF EVAPORATION   ***
!   ***                          OF RAIN                             ***
!   ***  COEFFS IS A COEFFICIENT GOVERNING THE RATE OF EVAPORATION   ***
!   ***                          OF SNOW                             ***
!   ***     CU IS THE COEFFICIENT GOVERNING CONVECTIVE MOMENTUM      ***
!   ***                         TRANSPORT                            ***
!   ***    DTMAX IS THE MAXIMUM NEGATIVE TEMPERATURE PERTURBATION    ***
!   ***        A LIFTED PARCEL IS ALLOWED TO HAVE BELOW ITS LFC      ***
!   ***    ALPHA AND DAMP ARE PARAMETERS THAT CONTROL THE RATE OF    ***
!   ***                 APPROACH TO QUASI-EQUILIBRIUM                ***
!   ***   (THEIR STANDARD VALUES ARE  0.20 AND 0.1, RESPECTIVELY)    ***
!   ***                   (DAMP MUST BE LESS THAN 1)                 ***

      integer, intent(in):: nd

! noff: integer limit for convection (nd-noff)
! minorig: First level of convection

      noff = 2
      minorig = 2

      nl=nd-noff
      nlp=nl+1
      nlm=nl-1

      elcrit=0.0011
      tlcrit=-55.0
      entp=1.5
      sigs=0.12
      sigd=0.05
      omtrain=50.0
      omtsnow=5.5
      coeffr=1.0
      coeffs=0.8
      dtmax=0.9
!
      cu=0.70
!
      betad=10.0
!
      damp=0.1
      alpha=0.2
!
      delta=0.01  ! cld
!
      return
      end
