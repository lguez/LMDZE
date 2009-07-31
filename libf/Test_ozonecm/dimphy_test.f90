module dimphy

  use dimens_m, only: jjm

  implicit none

  private jjm

  INTEGER, PARAMETER:: klon = jjm + 1

end module dimphy
