module dimens_m

  ! Model dimensions.

  implicit none

  INTEGER, PARAMETER:: iim = 16 ! number of longitudes
  INTEGER, PARAMETER:: jjm = 12 ! number of latitudes

  INTEGER, PARAMETER:: llm = 11
  ! number of vertical layers, should be >= 7, because of nl in cv30_param_m

  integer, parameter:: nqmx = 5 ! maximum number of tracers

end module dimens_m
