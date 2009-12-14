module parafilt

  ! From filtrez/parafilt.h, v 1.1.1.1 2004/05/19 12:53:09

  ! La premiere fois que  le Gcm  rentrera  dans le Filtre ,
  ! il indiquera  les bonnes valeurs  de  nfilun , nflius, nfilvn  et 
  ! nfilvs  a  mettre .  Il suffira alors de changer ces valeurs dans
  ! Parameter  ci-dessus  et de relancer  le  run.

  implicit none

  INTEGER, PARAMETER:: nfilun=6, nfilus=5, nfilvn=5, nfilvs=5

end module parafilt
