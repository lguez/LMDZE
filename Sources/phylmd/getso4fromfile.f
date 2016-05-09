module getso4fromfile_m

  implicit none

contains

  SUBROUTINE getso4fromfile(cyr, so4)

    ! Routine for reading SO4 data from files

    USE netcdf, ONLY: nf90_nowrite
    USE netcdf95, ONLY: nf95_close, nf95_get_var, nf95_inq_varid, &
         nf95_inquire_dimension, nf95_inquire_variable, nf95_open
    use nr_util, only: assert

    CHARACTER(len=*), intent(in):: cyr
    real, intent(out):: so4(:, :, :, :) ! (iim, jjm + 1, klev, 12)

    ! Local:
    CHARACTER(len=15) fname
    INTEGER NCID, VARID, nclen(4), i
    integer, pointer:: dimids(:) ! (4)

    !---------------------------------------------------------------------

    fname = 'aerosols' // cyr // '.nc'
    print *, 'Reading ', fname
    call NF95_OPEN(fname, NF90_NOWRITE, NCID)
    call NF95_INQ_VARID(NCID, "SO4", VARID)
    call nf95_inquire_variable(ncid, varid, dimids=dimids)

    do i = 1, 4
       call nf95_inquire_dimension(ncid, dimids(i), nclen=nclen(i))
    end do

    deallocate(dimids) ! pointer
    call assert(nclen == shape(so4), "getso4fromfile")
    call NF95_GET_VAR(NCID, VARID, so4)
    call NF95_CLOSE(NCID)

  END SUBROUTINE getso4fromfile

end module getso4fromfile_m
