module read_sst_m

  implicit none

contains

  SUBROUTINE read_sst(julien, knindex, lmt_sst)

    ! From interfoce_lim

    use conf_gcm_m, only: lmt_pas
    USE dimphy, ONLY: klon
    USE netcdf, ONLY: nf90_nowrite
    use netcdf95, only: NF95_CLOSE, nf95_get_var, NF95_INQ_VARID, nf95_open
    use nr_util, only: assert
    use time_phylmdz, only: itap

    integer, intent(IN):: julien ! jour a lire dans l'annee

    integer, intent(in):: knindex(:) ! (knon)
    ! index des points de la surface a traiter

    real, intent(out):: lmt_sst(:) ! (knon)
    ! SST lues dans le fichier de conditions aux limites

    ! Local:

    ! Champ lu dans le fichier de conditions aux limites :
    real, save:: sst_lu(klon)

    integer ncid, varid ! pour NetCDF

    ! --------------------------------------------------

    call assert(size(knindex) == size(lmt_sst), "read_sst knon")

    ! Tester d'abord si c'est le moment de lire le fichier
    if (mod(itap - 1, lmt_pas) == 0) then
       call NF95_OPEN ('limit.nc', NF90_NOWRITE, ncid)

       call NF95_INQ_VARID(ncid, 'SST', varid)
       call NF95_GET_VAR(ncid, varid, sst_lu, start = (/1, julien/))

       call NF95_CLOSE(ncid)
    endif

    lmt_sst = sst_lu(knindex)

  END SUBROUTINE read_sst

end module read_sst_m
