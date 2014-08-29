module cv_param

  ! From LMDZ4/libf/phylmd/cvparam.h,v 1.1.1.1 2004/05/19 12:53:08
  ! and LMDZ4/libf/phylmd/cv_routines.F, version 1.1.1.1, 2004/05/19 12:53:08

  ! Parameters for convectL: includes microphysical parameters, and
  ! parameters that control the rate of approach to quasi-equilibrium.

  use dimphy, only: klev

  implicit none

  integer, parameter:: noff= 2, minorig= 2, nl=klev-noff, nlp=nl+1, nlm=nl-1
  ! noff: integer limit for convection (klev-noff)
  ! minorig: First level of convection
  real, parameter:: elcrit=0.0011, tlcrit=-55.
  ! ELCRIT IS THE AUTOCONVERSION THERSHOLD WATER CONTENT (gm/gm) 
  ! TLCRIT IS CRITICAL TEMPERATURE BELOW WHICH THE AUTO- 
  ! CONVERSION THRESHOLD IS ASSUMED TO BE ZERO 
  ! (THE AUTOCONVERSION THRESHOLD VARIES LINEARLY 
  ! BETWEEN 0 C AND TLCRIT) 
  real, parameter:: entp=1.5
  ! ENTP IS THE COEFFICIENT OF MIXING IN THE ENTRAINMENT 
  ! FORMULATION 
  real, parameter:: sigs=0.12, sigd=0.05
  ! SIGD IS THE FRACTIONAL AREA COVERED BY UNSATURATED DNDRAFT 
  ! SIGS IS THE FRACTION OF PRECIPITATION FALLING OUTSIDE 
  ! OF CLOUD 
  real, parameter:: omtrain=50., omtsnow=5.5, coeffr=1., coeffs=0.8
  ! OMTRAIN IS THE ASSUMED FALL SPEED (P/s) OF RAIN 
  ! OMTSNOW IS THE ASSUMED FALL SPEED (P/s) OF SNOW 
  ! COEFFR IS A COEFFICIENT GOVERNING THE RATE OF EVAPORATION 
  ! OF RAIN 
  ! COEFFS IS A COEFFICIENT GOVERNING THE RATE OF EVAPORATION 
  ! OF SNOW 
  real, parameter:: dtmax=0.9
  ! DTMAX IS THE MAXIMUM NEGATIVE TEMPERATURE PERTURBATION 
  ! A LIFTED PARCEL IS ALLOWED TO HAVE BELOW ITS LFC 
  real, parameter:: cu=0.70
  ! CU IS THE COEFFICIENT GOVERNING CONVECTIVE MOMENTUM 
  ! TRANSPORT 
  real, parameter:: betad=10.
  real, parameter:: alpha=0.2, damp=0.1
  ! ALPHA AND DAMP ARE PARAMETERS THAT CONTROL THE RATE OF 
  ! APPROACH TO QUASI-EQUILIBRIUM 
  ! (THEIR STANDARD VALUES ARE 0.20 AND 0.1, RESPECTIVELY) 
  ! (DAMP MUST BE LESS THAN 1)
  real, parameter:: delta=0.01

end module cv_param
