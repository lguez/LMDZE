module ener

  ! From ener.h, v 1.1.1.1 2004/05/19 12:53:07 lmdzadmin

  use dimens_m, only: llm
  
  implicit none

  private llm

  REAL ang0, etot0, ptot0, ztot0, stot0, ang, etot, ptot, ztot, stot, rmsdpdt
  real rmsv, gtot(llm-1)

  save

end module ener
