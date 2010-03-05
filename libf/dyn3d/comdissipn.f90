module comdissipn

  ! From dyn3d/comdissipn.h,v 1.1.1.1 2004/05/19 12:53:06
  ! Les parametres de ce common proviennent des calculs effectues dans 
  ! Inidissip.

  use dimens_m, only: llm

  implicit none
  private llm

  real tetaudiv(llm),tetaurot(llm),tetah(llm)   
  real cdivu,      crot,         cdivh

end module comdissipn
