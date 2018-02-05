module indicesol

  implicit none

  INTEGER, parameter:: nbsrf=4

  INTEGER, parameter:: is_ter = 1 ! land
  INTEGER, parameter:: is_lic = 2 ! land ice
  INTEGER, parameter:: is_oce = 3 ! ocean
  INTEGER, parameter:: is_sic = 4 ! sea ice

  REAL, parameter:: epsfra = 1E-5
  CHARACTER(len=*), parameter:: clnsurf(nbsrf) = (/'ter', 'lic', 'oce', 'sic'/)

end module indicesol
