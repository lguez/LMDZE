module dimphy

  implicit none

  INTEGER, protected, save:: klon
  ! (number of distinct scalar grid points)
  ! (Les points de la physique sont les points scalaires de la dynamique.
  ! Numerotation :
  ! 1 pour le pole nord
  ! (jjm - 1) * iim pour l'interieur du domaine
  ! klon pour le pole sud)

  INTEGER, protected, save:: KLEV

contains

  subroutine init_dimphy

    use dimensions, only: iim, jjm, llm

    !--------------------------------------------------------

    klon = iim * (jjm - 1) + 2
    KLEV = llm

  end subroutine init_dimphy

end module dimphy
