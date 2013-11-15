module serre

  implicit none

  REAL:: clon = 0. ! longitude of the center of the zoom, in degrees
  real:: clat = 0. ! latitude of the center of the zoom, in degrees

  real, save:: transx, transy
  real, save:: alphax, alphay ! anciennes formulations des grossissements
  real, save:: pxo, pyo

  real:: grossismx = 1. ! facteur de grossissement du zoom, selon la longitude
  real:: grossismy = 1. ! facteur de grossissement du zoom, selon la latitude

  real:: dzoomx = 0.
  ! extension en longitude de la zone du zoom (fraction de la zone totale)

  real:: dzoomy = 0.
  ! extension en latitude de la zone du zoom (fraction de la zone totale)

  real:: taux = 3. ! raideur du zoom en X
  real:: tauy = 3. ! raideur du zoom en Y

end module serre
