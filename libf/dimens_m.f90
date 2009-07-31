module dimens_m

  ! Model dimensions.

  implicit none

  INTEGER, PARAMETER:: iim = 96 ! number of longitudes
  INTEGER, PARAMETER:: jjm = 72 ! number of latitudes
  INTEGER, PARAMETER:: llm = 50 ! number of vertical layers
  INTEGER, PARAMETER:: ndm = 1
  integer, parameter:: nqmx = 5 ! maximum number of tracers

end module dimens_m
