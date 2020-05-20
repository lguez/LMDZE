module dimensions

  ! Model dimensions.

  implicit none

  INTEGER, protected:: iim = 16 ! number of longitudes
  INTEGER, protected:: jjm = 12 ! number of latitudes

  INTEGER, protected:: llm = 11
  ! number of vertical layers, should be >= 7, because of nl in cv30_param_m

  integer, parameter:: nqmx = 5 ! maximum number of tracers

end module dimensions
