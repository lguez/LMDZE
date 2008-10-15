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
    ! Passage d'un champ de la grille physique à la grille dynamique

    use dimphy, only: klon

    REAL, intent(in):: pfi(:)
    real gr_fi_dyn(iim + 1, jjm + 1)

    ! Variable local to the procedure:
    real field(iim + 1, jjm + 1)

    !-----------------------------------------------------------------------

    if (size(pfi) /= klon) stop "gr_fi_dyn"

    ! Traitement des pôles :
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

  function gr_phy_write_2d(pfi)

    ! From phylmd/physiq.F, version 1.22 2006/02/20 09:38:28
    ! Transforme une variable de la grille physique à la grille d'écriture.
    ! The grid for output files does not duplicate the first longitude
    ! in the last longitude.

    use dimphy, only: klon

    REAL, intent(in):: pfi(:)
    real gr_phy_write_2d(iim, jjm + 1)

    ! Variable local to the procedure:
    real field(iim, jjm + 1)

    !-----------------------------------------------------------------------

    if (size(pfi) /= klon) stop "gr_phy_write_2d"

    ! Traitement des pôles :
    field(2:, 1) = pfi(1)
    field(2:, jjm + 1) = pfi(klon)

    gr_phy_write_2d = unpack(pfi, dyn_phy(:iim, :), field)

  END function gr_phy_write_2d

  !***************************************************

  function gr_phy_write_3d(pfi)

    ! Transforme une variable de la grille physique à la grille d'écriture.
    ! The grid for output files does not duplicate the first longitude
    ! in the last longitude.
    ! Input array has rank 2. Horizontal index is in the first dimension.

    use dimphy, only: klon
    use numer_rec, only: assert

    REAL, intent(in):: pfi(:, :)
    real gr_phy_write_3d(iim, jjm + 1, size(pfi, 2))

    ! Variable local to the procedure:
    integer l

    !-----------------------------------------------------------------------

    call assert(size(pfi, 1) == klon, "gr_phy_write_3d")

    do l = 1, size(pfi, 2)
       gr_phy_write_3d(:, :, l) = gr_phy_write_2d(pfi(:, l))
    end do

  END function gr_phy_write_3d

end module grid_change
