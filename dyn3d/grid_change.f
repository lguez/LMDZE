module grid_change

  use dimens_m, only: iim, jjm

  IMPLICIT NONE

  logical, save:: dyn_phy(iim + 1, jjm + 1)
  ! (mask for distinct points in the scalar grid and the "u" grid,
  ! first index is for longitude, second index is for latitude)

  private iim, jjm

contains

  subroutine init_dyn_phy

    ! Construct the mask:
    dyn_phy = .true.
    dyn_phy(2:, 1) = .false.
    dyn_phy(2:, jjm + 1) = .false.
    dyn_phy(iim + 1, 2:jjm) = .false.
    ! Note that "count(dyn_phy)" equals "klon"

  end subroutine init_dyn_phy

  !********************************************

  function gr_fi_dyn(pfi)

    ! From gr_fi_dyn.F, version 1.1.1.1 2004/05/19 12:53:05
    ! Passage d'un champ de la grille physique \`a la grille dynamique

    use dimphy, only: klon

    REAL, intent(in):: pfi(:) ! (klon)
    real gr_fi_dyn(iim + 1, jjm + 1)

    ! Variable local to the procedure:
    real field(iim + 1, jjm + 1)

    !-----------------------------------------------------------------------

    if (size(pfi) /= klon) stop "gr_fi_dyn"

    ! Traitement des p\^oles :
    field(2:, 1) = pfi(1)
    field(2:, jjm + 1) = pfi(klon)
    ! (We leave undefined elements in "field")

    gr_fi_dyn = unpack(pfi, dyn_phy, field)
    ! Undefined elements at last longitude in "gr_fi_dyn" come from
    ! undefined elements in "field".
    ! Overwrite them now, knowing that last longitude equals first longitude:
    gr_fi_dyn(iim + 1, 2:jjm) = gr_fi_dyn(1, 2:jjm)

  END function gr_fi_dyn

  !********************************************

  function gr_dyn_phy(v)

    ! Passage d'un champ 3D de la grille dynamique \`a la grille physique

    use dimphy, only: klon

    REAL, intent(in):: v(:, :, :) ! (iim + 1, jjm + 1, :)
    real gr_dyn_phy(klon, size(v, 3))

    ! Local:
    integer k

    !-----------------------------------------------------------------------

    forall (k = 1:size(v, 3)) gr_dyn_phy(:, k) = pack(v(:, :, k), dyn_phy)

  END function gr_dyn_phy

end module grid_change
