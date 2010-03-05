module coefils

  ! From filtrez/coefils.h,v 1.1.1.1 2004/05/19 12:53:09

  use dimens_m, only: iim, jjm

  implicit none

  private iim, jjm

  real sddu(iim),sddv(iim)
  real unsddu(iim),unsddv(iim),coefilu(iim,jjm),coefilv(iim,jjm)
  integer modfrstu(jjm),modfrstv(jjm)
  real eignfnu(iim,iim),eignfnv(iim,iim)
  real coefilu2(iim,jjm),coefilv2(iim,jjm)
  INTEGER jfiltnu,jfiltsu,jfiltnv,jfiltsv

end module coefils
