module comgeomphy

  ! From phylmd/comgeomphy.h,v 1.1.1.1 2004/05/19 12:53:09

  ! Passage de la geometrie de la dynamique a la physique

  use dimphy, only: klon

  implicit none

  real airephy(klon)
  real cuphy(klon), cvphy(klon) ! "cu" and "cv" values on the scalar grid
  REAL rlatd(klon), rlond(klon)

  save

  private klon

end module comgeomphy
