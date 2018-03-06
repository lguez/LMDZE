module dimphy

  use dimens_m, only: iim, jjm, llm

  implicit none

  private iim, jjm, llm

  INTEGER, PARAMETER:: klon = iim * (jjm - 1) + 2
  ! (number of distinct scalar grid points)
  ! (Les points de la physique sont les points scalaires de la dynamique.
  ! Numerotation :
  ! 1 pour le pole nord
  ! (jjm - 1) * iim pour l'interieur du domaine
  ! klon pour le pole sud)

  INTEGER, PARAMETER:: KLEV = llm
  REAL, save:: zmasq(KLON) ! fraction of land

end module dimphy