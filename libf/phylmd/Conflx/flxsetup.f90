      SUBROUTINE flxsetup
            use yoecumf
      IMPLICIT none
!
!     THIS ROUTINE DEFINES DISPOSABLE PARAMETERS FOR MASSFLUX SCHEME
!
!
      ENTRPEN=1.0E-4  ! ENTRAINMENT RATE FOR PENETRATIVE CONVECTION
      ENTRSCV=3.0E-4  ! ENTRAINMENT RATE FOR SHALLOW CONVECTION
      ENTRMID=1.0E-4  ! ENTRAINMENT RATE FOR MIDLEVEL CONVECTION
      ENTRDD =2.0E-4  ! ENTRAINMENT RATE FOR DOWNDRAFTS
      CMFCTOP=0.33  ! RELATIVE CLOUD MASSFLUX AT LEVEL ABOVE NONBUO LEVEL
      CMFCMAX=1.0  ! MAXIMUM MASSFLUX VALUE ALLOWED FOR UPDRAFTS ETC
      CMFCMIN=1.E-10  ! MINIMUM MASSFLUX VALUE (FOR SAFETY)
      CMFDEPS=0.3  ! FRACTIONAL MASSFLUX FOR DOWNDRAFTS AT LFS
      CPRCON =2.0E-4  ! CONVERSION FROM CLOUD WATER TO RAIN
      RHCDD=1.  ! RELATIVE SATURATION IN DOWNDRAFRS (NO LONGER USED)
!                 (FORMULATION IMPLIES SATURATION)
      LMFPEN = .TRUE.
      LMFSCV = .TRUE.
      LMFMID = .TRUE.
      LMFDD = .TRUE.
      LMFDUDV = .TRUE.
!
      RETURN
      END