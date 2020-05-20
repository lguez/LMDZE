module limit_read_sst_m

  implicit none

contains

  SUBROUTINE limit_read_sst(julien, knindex, tsurf)

    ! From interfoce_lim

    ! Libraries:
    USE netcdf, ONLY: nf90_nowrite
    use netcdf95, only: NF95_CLOSE, nf95_get_var, NF95_INQ_VARID, nf95_open
    use nr_util, only: assert

    use conf_gcm_m, only: lmt_pas
    USE dimphy, ONLY: klon
    use time_phylmdz, only: itap

    integer, intent(IN):: julien ! jour a lire dans l'annee

    integer, intent(in):: knindex(:) ! (knon)
    ! index des points de la surface a traiter

    real, intent(out):: tsurf(:) ! (knon)
    ! SST lues dans le fichier de conditions aux limites

    ! Local:

    ! Champ lu dans le fichier de conditions aux limites :
    real, save, allocatable:: sst_lu(:) ! (klon)

    integer ncid, varid ! pour NetCDF

    ! --------------------------------------------------

    call assert(size(knindex) == size(tsurf), "limit_read_sst knon")

    ! Tester d'abord si c'est le moment de lire le fichier
    if (mod(itap - 1, lmt_pas) == 0) then
       if (itap == 1) allocate(sst_lu(klon))
       call NF95_OPEN ('limit.nc', NF90_NOWRITE, ncid)

       call NF95_INQ_VARID(ncid, 'SST', varid)
       call NF95_GET_VAR(ncid, varid, sst_lu, start = (/1, julien/))

       call NF95_CLOSE(ncid)
    endif

    tsurf = sst_lu(knindex)

  END SUBROUTINE limit_read_sst

end module limit_read_sst_m
