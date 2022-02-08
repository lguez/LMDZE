module gr_phy_write_m

  use dimensions, only: iim, jjm
  use dimphy, only: klon

  IMPLICIT none

  interface gr_phy_write
     module procedure gr_phy_write_2d, gr_phy_write_3d
  end interface gr_phy_write

  private
  public gr_phy_write

contains

  function gr_phy_write_2d(pfi)

    ! From phylmd/physiq.F, version 1.22 2006/02/20 09:38:28
    ! Transforme une variable de la grille physique \`a la grille d'\'ecriture.
    ! The grid for output files does not duplicate the first longitude
    ! in the last longitude.

    use grid_change, only: dyn_phy

    REAL, intent(in):: pfi(:) ! (klon)
    real gr_phy_write_2d(iim, jjm + 1)

    ! Variable local to the procedure:
    real field(iim, jjm + 1)

    !-----------------------------------------------------------------------

    if (size(pfi) /= klon) stop "gr_phy_write_2d"

    ! Traitement des p\^oles :
    field(2:, 1) = pfi(1)
    field(2:, jjm + 1) = pfi(klon)

    gr_phy_write_2d = unpack(pfi, dyn_phy(:iim, :), field)

  END function gr_phy_write_2d

  !********************************************

  function gr_phy_write_3d(fi)

    ! From phylmd/physiq.F, version 1.22 2006/02/20 09:38:28

    ! Transforme une variable tri-dimensionnelle de la grille physique
    ! \`a la grille d'\'ecriture. The grid for output files does not
    ! duplicate the first longitude in the last longitude. Input array
    ! has rank 2. Horizontal index is in the first dimension.

    use jumble, only: assert

    REAL, intent(in):: fi(:, :) ! (klon, :)
    real gr_phy_write_3d(iim, jjm + 1, size(fi, 2))

    ! Local:
    INTEGER l

    !---------------

    call assert(size(fi, 1) == klon, "gr_phy_write_3d")

    DO l = 1, size(fi, 2)
       gr_phy_write_3d(:, :, l) = gr_phy_write_2d(fi(:, l))
    END DO

  END function gr_phy_write_3d

end module gr_phy_write_m
