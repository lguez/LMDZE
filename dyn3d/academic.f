module academic

  ! From dyn3d/academic.h,v 1.1.1.1 2004/05/19 12:53:05

  use dimens_m, only: iim, jjm, llm

  implicit none

  private iim, jjm, llm

  real tetarappel((iim + 1) * (jjm + 1), llm), taurappel

end module academic
