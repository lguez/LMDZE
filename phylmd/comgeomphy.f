module comgeomphy

  ! From phylmd/comgeomphy.h,v 1.1.1.1 2004/05/19 12:53:09

  ! Passage de la geometrie de la dynamique a la physique

  use dimphy, only: klon

  implicit none

  real, save:: airephy(klon)

  private klon

end module comgeomphy
