module dimensions

  ! Model dimensions.

  implicit none

  INTEGER, PARAMETER:: iim = CPP_IIM ! number of longitudes
  INTEGER, PARAMETER:: jjm = CPP_JJM ! number of latitudes

  INTEGER, PARAMETER:: llm = CPP_LLM
  ! number of vertical layers, should be >= 7, because of nl in cv30_param_m

  integer, parameter:: nqmx = 5 ! maximum number of tracers

end module dimensions
