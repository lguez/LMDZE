module raddim

  ! From phylmd/raddim.*.h, v 1.1.1.1 2004/05/19 12:53:08

  use dimens_m, only: llm
  use dimphy, only: klon

  implicit none

  INTEGER, PARAMETER:: kdlon=klon, kflev=llm

  private llm, klon

end module raddim
