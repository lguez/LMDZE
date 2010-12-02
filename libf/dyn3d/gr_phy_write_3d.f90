module gr_phy_write_3d_m

  implicit none

contains

  function gr_phy_write_3d(pfi)

    ! Transforme une variable de la grille physique à la grille d'écriture.
    ! The grid for output files does not duplicate the first longitude
    ! in the last longitude.
    ! Input array has rank 2. Horizontal index is in the first dimension.

    use dimens_m, only: iim, jjm
    use dimphy, only: klon
    use nr_util, only: assert
    use grid_change, only: gr_phy_write_2d

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

end module gr_phy_write_3d_m
