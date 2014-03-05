module yoecumf

  ! From phylmd/YOECUMF.h, version 1.1.1.1 2004/05/19 12:53:07

  ! M. Tiedtke, ECMWF, 18th Jan. 1989
  ! Parameters for cumulus massflux scheme

  implicit none

  logical lmfpen ! penetrative convection switched on
  logical lmfscv ! shallow convection switched on
  logical lmfmid ! midlevel convection switched on
  logical lmfdd ! cumulus downdraft switched on
  logical lmfdudv ! cumulus friction switched on
  real entrpen ! entrainment rate for penetrative convection
  real entrscv ! entrainment rate for shallow convection
  real entrmid ! entrainment rate for midlevel convection
  real entrdd ! entrainment rate for cumulus downdrafts
  real cmfctop ! relative cloud massflux at level above nonbuoyanc level
  real cmfcmax ! maximum massflux value allowed for updrafts etc
  real cmfcmin ! minimum massflux value (for safety)
  real cmfdeps ! fractional massflux for downdrafts at lfs

  real rhcdd
  ! relative saturation in downdrafts (no longer used) (formulation
  ! implies saturation)

  real cprcon ! coefficients for determining conversion from cloud water to rain

contains

  SUBROUTINE flxsetup

    ! This routine defines disposable parameters for massflux scheme.

    !--------------------------------------------------------

    ENTRPEN = 1E-4
    ENTRSCV = 3E-4
    ENTRMID = 1E-4
    ENTRDD = 2E-4 
    CMFCTOP = 0.33
    CMFCMAX = 1.
    CMFCMIN = 1E-10
    CMFDEPS = 0.3
    CPRCON = 2E-4
    RHCDD = 1.
    LMFPEN = .TRUE.
    LMFSCV = .TRUE.
    LMFMID = .TRUE.
    LMFDD = .TRUE.
    LMFDUDV = .TRUE.

  END SUBROUTINE flxsetup

end module yoecumf
