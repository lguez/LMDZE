module interfoce_lim_m

  implicit none

contains

  SUBROUTINE interfoce_lim(itime, dtime, jour, knindex, debut, lmt_sst, &
       pctsrf_new)

    ! lecture conditions limites
    ! Cette routine sert d'interface entre le modèle atmosphérique et
    ! un fichier de conditions aux limites.

    ! Laurent FAIRHEAD, February 2000

    use abort_gcm_m, only: abort_gcm
    USE dimphy, ONLY: klon
    USE indicesol, ONLY: is_lic, is_oce, is_sic, is_ter, nbsrf
    USE netcdf, ONLY: nf90_noerr, nf90_nowrite
    use netcdf95, only: NF95_CLOSE, nf95_get_var, NF95_INQ_VARID, nf95_open
    use nr_util, only: assert_eq

    integer, intent(IN):: itime ! numero du pas de temps courant
    real, intent(IN):: dtime ! pas de temps de la physique (en s)
    integer, intent(IN):: jour ! jour a lire dans l'annee

    integer, intent(in):: knindex(:) ! (knon)
    ! index des points de la surface a traiter

    logical, intent(IN):: debut ! 1er appel a la physique (initialisation)

    real, intent(out):: lmt_sst(:) ! (knon)
    ! SST lues dans le fichier de conditions aux limites

    real, intent(out):: pctsrf_new(:, :) ! (klon, nbsrf)
    ! sous-maille fractionnelle

    ! Local:

    integer knon ! nombre de points dans le domaine a traiter

    INTEGER, save:: lmt_pas ! frequence de lecture des conditions limites
    ! (en pas de physique)

    logical, save:: deja_lu
    ! pour indiquer que le jour à lire a déjà été lu pour une surface
    ! précédente

    integer, save:: jour_lu

    ! Champs lus dans le fichier de conditions aux limites :
    real, allocatable, save:: sst_lu(:)
    real, allocatable, save:: pct_tmp(:, :)

    integer ncid, varid ! pour NetCDF

    ! --------------------------------------------------

    knon = assert_eq(size(knindex), size(lmt_sst), "interfoce_lim knon")

    if (debut .and. .not. allocated(sst_lu)) then
       lmt_pas = nint(86400. / dtime) ! pour une lecture une fois par jour
       jour_lu = jour - 1
       allocate(sst_lu(klon))
       allocate(pct_tmp(klon, nbsrf))
    endif

    if ((jour - jour_lu) /= 0) deja_lu = .false.

    ! Tester d'abord si c'est le moment de lire le fichier
    if (mod(itime - 1, lmt_pas) == 0 .and. .not. deja_lu) then
       call NF95_OPEN ('limit.nc', NF90_NOWRITE, ncid)

       ! Fraction "ocean"
       call NF95_INQ_VARID(ncid, 'FOCE', varid)
       call NF95_GET_VAR(ncid, varid, pct_tmp(:, is_oce), start = (/1, jour/))

       ! Fraction "glace de mer"
       call NF95_INQ_VARID(ncid, 'FSIC', varid)
       call NF95_GET_VAR(ncid, varid, pct_tmp(:, is_sic), start = (/1, jour/))

       ! Fraction "terre"
       call NF95_INQ_VARID(ncid, 'FTER', varid)
       call NF95_GET_VAR(ncid, varid, pct_tmp(:, is_ter), start = (/1, jour/))

       ! Fraction "glacier terre"
       call NF95_INQ_VARID(ncid, 'FLIC', varid)
       call NF95_GET_VAR(ncid, varid, pct_tmp(:, is_lic), start = (/1, jour/))

       call NF95_INQ_VARID(ncid, 'SST', varid)
       call NF95_GET_VAR(ncid, varid, sst_lu, start = (/1, jour/))

       call NF95_CLOSE(ncid)
       deja_lu = .true.
       jour_lu = jour
    endif

    ! Recopie des variables dans les champs de sortie
    lmt_sst = sst_lu(knindex)
    pctsrf_new(:, is_oce) = pct_tmp(:, is_oce)
    pctsrf_new(:, is_sic) = pct_tmp(:, is_sic)

  END SUBROUTINE interfoce_lim

end module interfoce_lim_m
