module oasis_m

  ! From phylmd/oasis.h,v 1.1.1.1 2004/05/19 12:53:08

  !@  Contents : choice for the OASIS version: clim or pipe

  implicit none

  logical ok_oasis
  parameter(ok_oasis = .false.)

  CHARACTER(len=8) cchan
  PARAMETER ( cchan='CLIM' )

end module oasis_m
