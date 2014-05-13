module interfsur_lim_m

  implicit none

contains

  SUBROUTINE interfsur_lim(itime, dtime, jour, nisurf, knon, knindex, debut, &
       alb_new, z0_new)

    ! Cette routine sert d'interface entre le modèle atmosphérique et
    ! un fichier de conditions aux limites.

    ! L. Fairhead 02/2000

    use abort_gcm_m, only: abort_gcm
    USE dimphy, ONLY: klon
    use netcdf, only: NF90_NOWRITE
    use netcdf95, only: NF95_close, NF95_GET_VAR, NF95_INQ_VARID, NF95_OPEN

    integer, intent(IN):: itime ! numero du pas de temps courant
    real, intent(IN):: dtime ! pas de temps de la physique (en s)
    integer, intent(IN):: jour ! jour a lire dans l'annee

    integer, intent(IN):: nisurf
    ! index de la surface à traiter (1 = sol continental)

    integer, intent(IN):: knon ! nombre de points dans le domaine a traiter

    integer, intent(in):: knindex(:) ! (klon)
    ! index des points de la surface à traiter

    logical, intent(IN):: debut ! premier appel à la physique (initialisation)
    real, intent(out):: alb_new(:) ! (klon) albedo lu
    real, intent(out):: z0_new(:) ! (klon) longueur de rugosité lue

    ! Local:

    integer, save:: lmt_pas ! frequence de lecture des conditions limites
    ! (en pas de physique)

    logical, save:: deja_lu_sur
    ! jour à lire déjà lu pour une surface précédente

    integer, save:: jour_lu_sur

    ! Champs lus dans le fichier de conditions aux limites :
    real, allocatable, save:: alb_lu(:), rug_lu(:)

    integer ncid, varid

    !------------------------------------------------------------

    if (debut) then
       lmt_pas = nint(86400./dtime * 1.0) ! pour une lecture une fois par jour
       jour_lu_sur = jour - 1
       allocate(alb_lu(klon))
       allocate(rug_lu(klon))
    endif

    if (jour - jour_lu_sur /= 0) deja_lu_sur = .false.

    ! Tester d'abord si c'est le moment de lire le fichier
    if (mod(itime - 1, lmt_pas) == 0 .and. .not. deja_lu_sur) then
       call NF95_OPEN('limit.nc', NF90_NOWRITE, ncid)

       ! Lecture Albedo
       call NF95_INQ_VARID(ncid, 'ALB', varid)
       call NF95_GET_VAR(ncid, varid, alb_lu, start=(/1, jour/))

       ! Lecture rugosité
       call NF95_INQ_VARID(ncid, 'RUG', varid)
       call NF95_GET_VAR(ncid, varid, rug_lu, start=(/1, jour/))

       call NF95_CLOSE(ncid)
       deja_lu_sur = .true.
       jour_lu_sur = jour
    endif

    ! Recopie des variables dans les champs de sortie
    alb_new(:knon) = alb_lu(knindex(:knon))
    z0_new(:knon) = rug_lu(knindex(:knon))
    alb_new(knon + 1:) = 999999.
    z0_new(knon + 1:) = 999999.

  END SUBROUTINE interfsur_lim

end module interfsur_lim_m
