module serre

  implicit none

  REAL:: clon= 0. ! longitude en degres du centre du zoom

  real:: clat= 0. ! latitude en degres du centre du zoom

  real, save:: transx,transy
  real, save:: alphax, alphay ! anciennes formulations des grossissements
  real, save:: pxo,pyo

  real:: grossismx= 1.0 ! facteur de grossissement du zoom, selon la longitude

  real:: grossismy= 1.0 ! facteur de grossissement du zoom, selon la latitude

  real:: dzoomx= 0.0
  ! (extension en longitude de la zone du zoom (fraction de la zone totale))

  real:: dzoomy= 0.0
  ! extension en latitude de la zone du zoom 
  ! (fraction de la zone totale)

  real:: taux= 3.0 ! raideur du zoom en X

  real:: tauy= 3.0 ! raideur du zoom en Y

end module serre
