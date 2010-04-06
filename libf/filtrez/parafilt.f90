module parafilt

  ! From filtrez/parafilt.h, v 1.1.1.1 2004/05/19 12:53:09

  ! La premiere fois que  le Gcm  rentrera  dans le Filtre ,
  ! il indiquera  les bonnes valeurs  de  nfilun , nflius, nfilvn  et 
  ! nfilvs  a  mettre . Il suffira alors de changer ces valeurs dans
  ! Parameter  ci-dessus et de relancer  le  run.

  use dimens_m, only: iim

  implicit none

  private iim

  INTEGER, PARAMETER:: nfilun=3, nfilus=2, nfilvn=2, nfilvs=2

  real matriceun(iim,iim,nfilun), matriceus(iim,iim,nfilus)
  real matricevn(iim,iim,nfilvn), matricevs(iim,iim,nfilvs)
  real matrinvn(iim,iim,nfilun), matrinvs(iim,iim,nfilus)

end module parafilt
