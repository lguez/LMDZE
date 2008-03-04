module dimens_m

  ! Model dimensions.

  implicit none

  INTEGER, PARAMETER:: iim = 32 ! number of longitudes
  INTEGER, PARAMETER:: jjm = 24 ! number of latitudes
  INTEGER, PARAMETER:: llm = 9 ! number of vertical layers
  INTEGER, PARAMETER:: ndm = 1
  integer, parameter:: nqmx = 4 ! maximum number of tracers

end module dimens_m
