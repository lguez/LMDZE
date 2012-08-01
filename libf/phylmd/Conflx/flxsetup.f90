module flxsetup_m

  IMPLICIT none

contains

  SUBROUTINE flxsetup

    ! This routine defines disposable parameters for massflux scheme.

    USE yoecumf, ONLY: cmfcmax, cmfcmin, cmfctop, cmfdeps, cprcon, entrdd, &
         entrmid, entrpen, entrscv, lmfdd, lmfdudv, lmfmid, lmfpen, lmfscv, &
         rhcdd

    !--------------------------------------------------------

    ENTRPEN=1E-4 ! ENTRAINMENT RATE FOR PENETRATIVE CONVECTION
    ENTRSCV=3E-4 ! ENTRAINMENT RATE FOR SHALLOW CONVECTION
    ENTRMID=1E-4 ! ENTRAINMENT RATE FOR MIDLEVEL CONVECTION
    ENTRDD =2E-4 ! ENTRAINMENT RATE FOR DOWNDRAFTS
    CMFCTOP=0.33 ! RELATIVE CLOUD MASSFLUX AT LEVEL ABOVE NONBUO LEVEL
    CMFCMAX=1. ! MAXIMUM MASSFLUX VALUE ALLOWED FOR UPDRAFTS ETC
    CMFCMIN=1.E-10 ! MINIMUM MASSFLUX VALUE (FOR SAFETY)
    CMFDEPS=0.3 ! FRACTIONAL MASSFLUX FOR DOWNDRAFTS AT LFS
    CPRCON =2E-4 ! CONVERSION FROM CLOUD WATER TO RAIN

    RHCDD=1.
    ! RELATIVE SATURATION IN DOWNDRAFRS (NO LONGER USED) (FORMULATION
    ! IMPLIES SATURATION)

    LMFPEN = .TRUE.
    LMFSCV = .TRUE.
    LMFMID = .TRUE.
    LMFDD = .TRUE.
    LMFDUDV = .TRUE.

  END SUBROUTINE flxsetup

end module flxsetup_m
