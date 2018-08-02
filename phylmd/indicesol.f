module indicesol

  implicit none

  INTEGER, parameter:: nbsrf=4

  INTEGER, parameter:: is_ter = 1 ! land
  ! Note: the first surface type must be land. pbl_surface must call
  ! the land interface first because the ocean needs the run-off.

  INTEGER, parameter:: is_lic = 2 ! land ice

  INTEGER, parameter:: is_oce = 3 ! ocean
  ! Ocean must be before sea ice.
  
  INTEGER, parameter:: is_sic = 4 ! sea ice

  REAL, parameter:: epsfra = 1E-5
  CHARACTER(len=*), parameter:: clnsurf(nbsrf) = (/'ter', 'lic', 'oce', 'sic'/)

end module indicesol
