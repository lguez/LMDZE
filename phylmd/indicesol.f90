module indicesol

  implicit none

  INTEGER, parameter:: nbsrf=4
  INTEGER, parameter:: is_oce=3 ! ocean
  INTEGER, parameter:: is_sic = 4 ! sea ice
  INTEGER, parameter:: is_ter = 1 ! land
  INTEGER, parameter:: is_lic = 2 ! land ice
  REAL, parameter:: epsfra = 1E-05
  CHARACTER(len=*), parameter:: clnsurf(nbsrf) = (/'ter', 'lic', 'oce', 'sic'/)

end module indicesol
