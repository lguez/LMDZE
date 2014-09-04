module interfoce_lim_m

  implicit none

contains

  SUBROUTINE interfoce_lim(itime, dtime, jour, klon, knon, knindex, &
       debut, lmt_sst, pctsrf_new)

    ! Cette routine sert d'interface entre le modèle atmosphérique et
    ! un fichier de conditions aux limites.

    ! Laurent FAIRHEAD, February 2000

    use abort_gcm_m, only: abort_gcm
    USE indicesol, ONLY: is_lic, is_oce, is_sic, is_ter, nbsrf
    USE netcdf, ONLY: nf90_noerr, nf90_nowrite
    use netcdf95, only: NF95_CLOSE, nf95_get_var, NF95_INQ_VARID, nf95_open

    integer, intent(IN):: itime ! numero du pas de temps courant
    real, intent(IN):: dtime ! pas de temps de la physique (en s)
    integer, intent(IN):: jour ! jour a lire dans l'annee
    integer, intent(IN):: klon ! taille de la grille
    integer, intent(IN):: knon ! nombre de points dans le domaine a traiter

    integer, intent(in):: knindex(klon)
    ! index des points de la surface a traiter

    logical, intent(IN):: debut ! 1er appel a la physique (initialisation)

    real, intent(out), dimension(knon):: lmt_sst
    ! SST lues dans le fichier de CL

    real, intent(out), dimension(klon, nbsrf):: pctsrf_new
    ! sous-maille fractionnelle

    ! Local:

    integer ii

    INTEGER, save:: lmt_pas ! frequence de lecture des conditions limites
    ! (en pas de physique)

    logical, save:: deja_lu ! pour indiquer que le jour a lire a deja
    ! lu pour une surface precedente

    integer, save:: jour_lu

    ! Champs lus dans le fichier de conditions aux limites :
    real, allocatable, save, dimension(:):: sst_lu
    real, allocatable, save, dimension(:, :):: pct_tmp

    ! Pour Netcdf :
    integer nid, nvarid
    integer, dimension(2):: start, epais

    ! --------------------------------------------------

    if (debut .and. .not. allocated(sst_lu)) then
       lmt_pas = nint(86400./dtime * 1.0) ! pour une lecture une fois par jour
       jour_lu = jour - 1
       allocate(sst_lu(klon))
       allocate(pct_tmp(klon, nbsrf))
    endif

    if ((jour - jour_lu) /= 0) deja_lu = .false.

    ! Tester d'abord si c'est le moment de lire le fichier
    if (mod(itime - 1, lmt_pas) == 0 .and. .not. deja_lu) then
       call NF95_OPEN ('limit.nc', NF90_NOWRITE, nid)

       ! La tranche de donnees a lire:
       start(1) = 1
       start(2) = jour
       epais(1) = klon
       epais(2) = 1

       ! Fraction "ocean"
       call NF95_INQ_VARID(nid, 'FOCE', nvarid)
       call NF95_GET_VAR(nid, nvarid, pct_tmp(:, is_oce), start, epais)

       ! Fraction "glace de mer"
       call NF95_INQ_VARID(nid, 'FSIC', nvarid)
       call NF95_GET_VAR(nid, nvarid, pct_tmp(:, is_sic), start, epais)

       ! Fraction "terre"
       call NF95_INQ_VARID(nid, 'FTER', nvarid)
       call NF95_GET_VAR(nid, nvarid, pct_tmp(:, is_ter), start, epais)

       ! Fraction "glacier terre"
       call NF95_INQ_VARID(nid, 'FLIC', nvarid)
       call NF95_GET_VAR(nid, nvarid, pct_tmp(:, is_lic), start, epais)

       call NF95_INQ_VARID(nid, 'SST', nvarid)
       call NF95_GET_VAR(nid, nvarid, sst_lu, start, epais)

       call NF95_CLOSE(nid)
       deja_lu = .true.
       jour_lu = jour
    endif

    ! Recopie des variables dans les champs de sortie

    do ii = 1, knon
       lmt_sst(ii) = sst_lu(knindex(ii))
    enddo

    pctsrf_new(:, is_oce) = pct_tmp(:, is_oce)
    pctsrf_new(:, is_sic) = pct_tmp(:, is_sic)

  END SUBROUTINE interfoce_lim

end module interfoce_lim_m
