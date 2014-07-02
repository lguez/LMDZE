module gath2cpl_m

  ! From phylmd/interface_surf.F90, version 1.8, 2005/05/25 13:10:09

  implicit none

contains

  SUBROUTINE gath2cpl(champ_in, champ_out, klon, knon, iim, jjm, knindex)

    ! Cette routine ecrit un champ 'gathered' sur la grille 2D pour le passer
    ! au coupleur.

    integer, intent(in):: klon, knon, iim, jjm
    ! klon taille de la grille
    ! knon nombre de points dans le domaine a traiter
    ! iim, jjm dimension de la grille 2D

    real, intent(in):: champ_in(klon) ! champ sur la grille gathere 

    integer, intent(in):: knindex(klon)
    ! index des points de la surface a traiter

    real, intent(out):: champ_out(iim, jjm + 1) ! champ sur la grille 2D

    ! Local:
    integer :: i, ig, j
    real, dimension(klon) :: tamp

    !--------------------------------------------------------------------

    tamp = 0.
    do i = 1, knon
       ig = knindex(i)
       tamp(ig) = champ_in(i)
    enddo
    ig = 1
    champ_out(:, 1) = tamp(ig)
    do j = 2, jjm
       do i = 1, iim
          ig = ig + 1
          champ_out(i, j) = tamp(ig)
       enddo
    enddo
    ig = ig + 1
    champ_out(:, jjm + 1) = tamp(ig)

  END SUBROUTINE gath2cpl

end module gath2cpl_m
