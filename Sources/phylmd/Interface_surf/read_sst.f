module read_sst_m

  implicit none

contains

  SUBROUTINE read_sst(dtime, jour, knindex, debut, lmt_sst)

    ! From interfoce_lim

    USE dimphy, ONLY: klon
    USE netcdf, ONLY: nf90_nowrite
    use netcdf95, only: NF95_CLOSE, nf95_get_var, NF95_INQ_VARID, nf95_open
    use nr_util, only: assert
    use time_phylmdz, only: itap

    real, intent(IN):: dtime ! pas de temps de la physique (en s)
    integer, intent(IN):: jour ! jour a lire dans l'annee

    integer, intent(in):: knindex(:) ! (knon)
    ! index des points de la surface a traiter

    logical, intent(IN):: debut ! 1er appel a la physique (initialisation)

    real, intent(out):: lmt_sst(:) ! (knon)
    ! SST lues dans le fichier de conditions aux limites

    ! Local:

    INTEGER, save:: lmt_pas ! frequence de lecture des conditions limites
    ! (en pas de physique)

    logical, save:: deja_lu
    ! pour indiquer que le jour à lire a déjà été lu pour une surface
    ! précédente

    integer, save:: jour_lu

    ! Champ lu dans le fichier de conditions aux limites :
    real, allocatable, save:: sst_lu(:)

    integer ncid, varid ! pour NetCDF

    ! --------------------------------------------------

    call assert(size(knindex) == size(lmt_sst), "read_sst knon")

    if (debut .and. .not. allocated(sst_lu)) then
       lmt_pas = nint(86400. / dtime) ! pour une lecture une fois par jour
       jour_lu = jour - 1
       allocate(sst_lu(klon))
    endif

    if ((jour - jour_lu) /= 0) deja_lu = .false.

    ! Tester d'abord si c'est le moment de lire le fichier
    if (mod(itap - 1, lmt_pas) == 0 .and. .not. deja_lu) then
       call NF95_OPEN ('limit.nc', NF90_NOWRITE, ncid)

       call NF95_INQ_VARID(ncid, 'SST', varid)
       call NF95_GET_VAR(ncid, varid, sst_lu, start = (/1, jour/))

       call NF95_CLOSE(ncid)
       deja_lu = .true.
       jour_lu = jour
    endif

    ! Recopie dans le champ de sortie
    lmt_sst = sst_lu(knindex)

  END SUBROUTINE read_sst

end module read_sst_m
