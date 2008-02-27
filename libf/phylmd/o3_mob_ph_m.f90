module o3_Mob_ph_m

  implicit none

contains

  function o3_Mob_ph(ncid, name)

    ! This function reads a single Mobidic ozone parameter from a file and
    ! packs it on the "physics" grid.

    use dimens_m, only: iim, jjm, llm
    use dimphy, only: klon
    use netcdf95, only: nf95_inq_varid, nf90_get_var, handle_err
    use grid_change, only: dyn_phy

    integer, intent(in):: ncid ! NetCDF ID of the file
    character(len=*), intent(in):: name ! of the NetCDF variable

    real o3_Mob_ph(klon, llm, 12)
    ! (ozone parameter from Mobidic on the "physics" grid)
    ! (Third dimension is the number of the month in the year.
    ! "o3_Mob_ph(i, k, month)" is at longitude "xlon(i)", latitude
    ! "xlat(i)", middle of layer "k".)

    ! Variables local to the procedure:
    integer varid, ncerr
    integer k, month

    real o3_Mob_dyn(iim + 1, jjm + 1, llm, 12)
    ! (ozone parameter from Mobidic on the "dynamics" grid)
    ! Fourth dimension is the number of the month in the year.
    ! "o3_Mob_dyn(i, j, k, month)" is at longitude "rlonv(i)", latitude
    ! "rlatu(j)", middle of layer "k".)

    !--------------------------------------------

    call nf95_inq_varid(ncid, name, varid)
    ncerr = nf90_get_var(ncid, varid, o3_Mob_dyn)
    call handle_err("o3_Mob_ph nf90_get_var " // name, ncerr, ncid)

    ! Latitudes are in increasing order in the input file while
    ! "rlatu" is in decreasing order, so invert:
    o3_Mob_dyn = o3_Mob_dyn(:, jjm+1:1:-1, :, :)
    forall (k = 1:llm, month = 1:12) &
         o3_Mob_ph(:, k, month) = pack(o3_Mob_dyn(:, :, k, month), dyn_phy)

  end function o3_Mob_ph

end module o3_Mob_ph_m
