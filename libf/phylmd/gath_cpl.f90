module gath_cpl

  ! From phylmd/interface_surf.F90,v 1.8 2005/05/25 13:10:09

  implicit none

contains

  SUBROUTINE gath2cpl(champ_in, champ_out, klon, knon, iim, jjm, knindex)

    ! Cette routine ecrit un champ 'gathered' sur la grille 2D pour le passer
    ! au coupleur.
    !
    ! 
    ! input:         
    !   champ_in     champ sur la grille gathere        
    !   knon         nombre de points dans le domaine a traiter
    !   knindex      index des points de la surface a traiter
    !   klon         taille de la grille
    !   iim,jjm      dimension de la grille 2D
    !
    ! output:
    !   champ_out    champ sur la grille 2D
    !
    ! input
    integer                   :: klon, knon, iim, jjm
    real, dimension(klon)     :: champ_in
    integer, dimension(klon)  :: knindex
    ! output
    real, dimension(iim,jjm+1)  :: champ_out
    ! local
    integer                   :: i, ig, j
    real, dimension(klon)     :: tamp

    tamp = 0.
    do i = 1, knon
       ig = knindex(i)
       tamp(ig) = champ_in(i)
    enddo
    ig = 1
    champ_out(:,1) = tamp(ig)
    do j = 2, jjm
       do i = 1, iim
          ig = ig + 1
          champ_out(i,j) = tamp(ig)
       enddo
    enddo
    ig = ig + 1
    champ_out(:,jjm+1) = tamp(ig)

  END SUBROUTINE gath2cpl

  !*****************************************************

  SUBROUTINE cpl2gath(champ_in, champ_out, klon, knon, iim, jjm, knindex)

    ! Cette routine ecrit un champ 'gathered' sur la grille 2D pour le passer
    ! au coupleur.
    !
    ! 
    ! input:         
    !   champ_in     champ sur la grille gathere        
    !   knon         nombre de points dans le domaine a traiter
    !   knindex      index des points de la surface a traiter
    !   klon         taille de la grille
    !   iim,jjm      dimension de la grille 2D
    !
    ! output:
    !   champ_out    champ sur la grille 2D
    !
    ! input
    integer                   :: klon, knon, iim, jjm
    real, dimension(iim,jjm+1)     :: champ_in
    integer, dimension(klon)  :: knindex
    ! output
    real, dimension(klon)  :: champ_out
    ! local
    integer                   :: i, ig, j
    real, dimension(klon)     :: tamp

    ig = 1
    tamp(ig) = champ_in(1,1)
    do j = 2, jjm
       do i = 1, iim
          ig = ig + 1
          tamp(ig) = champ_in(i,j)
       enddo
    enddo
    ig = ig + 1
    tamp(ig) = champ_in(1,jjm+1)

    do i = 1, knon
       ig = knindex(i)
       champ_out(i) = tamp(ig)
    enddo

  END SUBROUTINE cpl2gath

end module gath_cpl
