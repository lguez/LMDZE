module regr_lat_time_coefoz_m

  ! This module is clean: no C preprocessor directive, no include line.
  ! Author: Lionel GUEZ

  implicit none

  private
  public regr_lat_time_coefoz

contains

  subroutine regr_lat_time_coefoz

    ! "regr_lat_time_coefoz" stands for "regrid latitude time
    ! coefficients ozone".

    ! This procedure reads from a NetCDF file parameters for ozone
    ! chemistry, regrids them in latitude and time, and writes the
    ! regridded fields to a new NetCDF file.

    ! The input fields depend on time, pressure level and
    ! latitude.
    ! We assume that the input fields are step functions
    ! of latitude.
    ! Regridding in latitude is made by averaging, with a cosine of
    ! latitude factor.
    ! The target LMDZ latitude grid is the "scalar" grid: "rlatu".
    ! The values of "rlatu" are taken to be the centers of intervals.
    ! Regridding in time is by linear interpolation.
    ! Monthly values are processed to get daily values, on the basis
    ! of a 360-day calendar.

    ! We assume that in the input file:
    ! -- the latitude is in degrees and strictly monotonic (as all
    ! NetCDF coordinate variables should be);
    ! -- time increases from January to December (even though we do
    ! not use values of the input time coordinate).

    use dimens_m, only: jjm
    use comgeom, only: rlatv
    use nr_util, only: pi
    use regr1_step_av_m, only: regr1_step_av
    use regr3_lint_m, only: regr3_lint
    use netcdf95, only: nf95_open, nf95_gw_var, nf95_close, &
         nf95_inq_varid, handle_err, nf95_put_var
    use netcdf, only: nf90_nowrite, nf90_get_var

    ! Variables local to the procedure:

    integer ncid_in, ncid_out ! NetCDF IDs for input and output files
    integer n_plev ! number of pressure levels in the input data
    integer n_lat! number of latitudes in the input data

    real, pointer:: latitude(:)
    ! (of input data, converted to rad, sorted in strictly increasing order)

    real, allocatable:: lat_in_edg(:)
    ! (edges of latitude intervals for input data, in rad, in strictly
    ! increasing order)

    real, pointer:: plev(:) ! pressure level of input data
    logical decr_lat ! decreasing latitude in the input file

    real, allocatable:: o3_par_in(:, :, :) ! (n_lat, n_plev, 12)
    ! (ozone parameter from the input file)
    ! ("o3_par_in(j, l, month)" is at latitude "latitude(j)" and pressure
    ! level "plev(l)". "month" is between 1 and 12.)

    real, allocatable:: v_regr_lat(:, :, :) ! (jjm + 1, n_plev, 0:13)
    ! (mean of a variable "v" over a latitude interval)
    ! (First dimension is latitude interval.
    ! The latitude interval for "v_regr_lat(j,:, :)" contains "rlatu(j)".
    ! If "j" is between 2 and "jjm" then the interval is:
    ! [rlatv(j), rlatv(j-1)]
    ! If "j" is 1 or "jjm + 1" then the interval is:
    ! [rlatv(1), pi / 2]
    ! or:
    ! [- pi / 2, rlatv(jjm)]
    ! respectively.
    ! "v_regr_lat(:, l, :)" is for pressure level "plev(l)".
    ! Last dimension is month number.)

    real, allocatable:: o3_par_out(:, :, :) ! (jjm + 1, n_plev, 360)
    ! (regridded ozone parameter)
    ! ("o3_par_out(j, l, day)" is at latitude "rlatu(j)", pressure
    ! level "plev(l)" and date "January 1st 0h" + "tmidday(day)", in a
    ! 360-day calendar.)

    integer j
    integer i_v ! index of ozone parameter
    integer, parameter:: n_o3_param = 8 ! number of ozone parameters

    character(len=11) name_in(n_o3_param)
    ! (name of NetCDF primary variable in the input file)

    character(len=9) name_out(n_o3_param)
    ! (name of NetCDF primary variable in the output file)

    integer varid_in(n_o3_param), varid_out(n_o3_param), varid_plev, varid_time
    integer ncerr, varid
    ! (for NetCDF)

    real, parameter:: tmidmonth(0:13) = (/(-15. + 30. * j, j = 0, 13)/)
    ! (time to middle of month, in days since January 1st 0h, in a
    ! 360-day calendar)
    ! (We add values -15 and 375 so that, for example, day 3 of the year is
    ! interpolated between the December and the January value.)

    real, parameter:: tmidday(360) = (/(j + 0.5, j = 0, 359)/)
    ! (time to middle of day, in days since January 1st 0h, in a
    ! 360-day calendar)

    !---------------------------------

    print *, "Call sequence information: regr_lat_time_coefoz"

    ! Names of ozone parameters:
    i_v = 0

    i_v = i_v + 1
    name_in(i_v) = "P_net"
    name_out(i_v) = "P_net_Mob"

    i_v = i_v + 1
    name_in(i_v) = "a2"
    name_out(i_v) = "a2"

    i_v = i_v + 1
    name_in(i_v) = "tro3"
    name_out(i_v) = "r_Mob"

    i_v = i_v + 1
    name_in(i_v) = "a4"
    name_out(i_v) = "a4"

    i_v = i_v + 1
    name_in(i_v) = "temperature"
    name_out(i_v) = "temp_Mob"

    i_v = i_v + 1
    name_in(i_v) = "a6"
    name_out(i_v) = "a6"

    i_v = i_v + 1
    name_in(i_v) = "Sigma"
    name_out(i_v) = "Sigma_Mob"

    i_v = i_v + 1
    name_in(i_v) = "R_Het"
    name_out(i_v) = "R_Het"

    call nf95_open("coefoz.nc", nf90_nowrite, ncid_in)

    ! Get coordinates from the input file:

    call nf95_inq_varid(ncid_in, "latitude", varid)
    call nf95_gw_var(ncid_in, varid, latitude)
    ! Convert from degrees to rad, because "rlatv" is in rad:
    latitude = latitude / 180. * pi
    n_lat = size(latitude)
    ! We need to supply the latitudes to "regr1_step_av" in
    ! increasing order, so invert order if necessary:
    decr_lat = latitude(1) > latitude(n_lat)
    if (decr_lat) latitude = latitude(n_lat:1:-1)

    ! Compute edges of latitude intervals:
    allocate(lat_in_edg(n_lat + 1))
    lat_in_edg(1) = - pi / 2
    forall (j = 2:n_lat) lat_in_edg(j) = (latitude(j - 1) + latitude(j)) / 2
    lat_in_edg(n_lat + 1) = pi / 2
    deallocate(latitude) ! pointer

    call nf95_inq_varid(ncid_in, "plev", varid)
    call nf95_gw_var(ncid_in, varid, plev)
    n_plev = size(plev)
    ! (We only need the pressure coordinate to copy it to the output file.)

    ! Get the IDs of ozone parameters in the input file:
    do i_v = 1, n_o3_param
       call nf95_inq_varid(ncid_in, trim(name_in(i_v)), varid_in(i_v))
    end do

    ! Create the output file and get the variable IDs:
    call prepare_out(ncid_in, varid_in, n_plev, name_out, ncid_out, &
         varid_out, varid_plev, varid_time)

    ! Write remaining coordinate variables:
    call nf95_put_var(ncid_out, varid_time, tmidday)
    call nf95_put_var(ncid_out, varid_plev, plev)

    deallocate(plev) ! pointer

    allocate(o3_par_in(n_lat, n_plev, 12))
    allocate(v_regr_lat(jjm + 1, n_plev, 0:13))
    allocate(o3_par_out(jjm + 1, n_plev, 360))

    do i_v = 1, n_o3_param
       ! Process ozone parameter "name_in(i_v)"

       ncerr = nf90_get_var(ncid_in, varid_in(i_v), o3_par_in)
       call handle_err("nf90_get_var", ncerr, ncid_in)

       if (decr_lat) o3_par_in = o3_par_in(n_lat:1:-1, :, :)

       ! Regrid in latitude:
       ! We average with respect to sine of latitude, which is
       ! equivalent to weighting by cosine of latitude:
       v_regr_lat(jjm+1:1:-1, :, 1:12) = regr1_step_av(o3_par_in, &
            xs=sin(lat_in_edg), xt=sin((/- pi / 2, rlatv(jjm:1:-1), pi / 2/)))
       ! (invert order of indices in "v_regr_lat" because "rlatu" is
       ! decreasing)

       ! Duplicate January and December values, in preparation of time
       ! interpolation:
       v_regr_lat(:, :, 0) = v_regr_lat(:, :, 12)
       v_regr_lat(:, :, 13) = v_regr_lat(:, :, 1)

       ! Regrid in time by linear interpolation:
       o3_par_out = regr3_lint(v_regr_lat, tmidmonth, tmidday)

       ! Write to file:
       call nf95_put_var(ncid_out, varid_out(i_v), &
            o3_par_out(jjm+1:1:-1, :, :))
       ! (The order of "rlatu" is inverted in the output file)
    end do

    call nf95_close(ncid_out)
    call nf95_close(ncid_in)

  end subroutine regr_lat_time_coefoz

  !********************************************

  subroutine prepare_out(ncid_in, varid_in, n_plev, name_out, ncid_out, &
       varid_out, varid_plev, varid_time)

    ! This subroutine creates the NetCDF output file, defines
    ! dimensions and variables, and writes one of the coordinate variables.

    use dimens_m, only: jjm
    use comgeom, only: rlatu
    use nr_util, only: assert_eq, pi

    use netcdf95, only: nf95_create, nf95_def_dim, nf95_def_var, &
         nf95_put_att, nf95_enddef, nf95_copy_att, nf95_put_var
    use netcdf, only: nf90_clobber, nf90_float, nf90_copy_att, nf90_global

    integer, intent(in):: ncid_in, varid_in(:), n_plev
    character(len=*), intent(in):: name_out(:) ! of NetCDF variables
    integer, intent(out):: ncid_out, varid_out(:), varid_plev, varid_time

    ! Variables local to the procedure:
    integer ncerr
    integer dimid_rlatu, dimid_plev, dimid_time
    integer varid_rlatu
    integer i, n_o3_param

    !---------------------------

    print *, "Call sequence information: prepare_out"
    n_o3_param = assert_eq(size(varid_in), size(varid_out), &
         size(name_out), "prepare_out") 

    call nf95_create("coefoz_LMDZ.nc", nf90_clobber, ncid_out)

    ! Dimensions:
    call nf95_def_dim(ncid_out, "time", 360, dimid_time)
    call nf95_def_dim(ncid_out, "plev", n_plev, dimid_plev)
    call nf95_def_dim(ncid_out, "rlatu", jjm + 1, dimid_rlatu)

    ! Define coordinate variables:

    call nf95_def_var(ncid_out, "time", nf90_float, dimid_time, varid_time)
    call nf95_put_att(ncid_out, varid_time, "units", "days since 2000-1-1")
    call nf95_put_att(ncid_out, varid_time, "calendar", "360_day")
    call nf95_put_att(ncid_out, varid_time, "standard_name", "time")

    call nf95_def_var(ncid_out, "plev", nf90_float, dimid_plev, varid_plev)
    call nf95_put_att(ncid_out, varid_plev, "units", "millibar")
    call nf95_put_att(ncid_out, varid_plev, "standard_name", "air_pressure")
    call nf95_put_att(ncid_out, varid_plev, "long_name", "air pressure")

    call nf95_def_var(ncid_out, "rlatu", nf90_float, dimid_rlatu, varid_rlatu)
    call nf95_put_att(ncid_out, varid_rlatu, "units", "degrees_north")
    call nf95_put_att(ncid_out, varid_rlatu, "standard_name", "latitude")

    ! Define primary variables:

    do i = 1, n_o3_param
       call nf95_def_var(ncid_out, name_out(i), nf90_float, &
            (/dimid_rlatu, dimid_plev, dimid_time/), varid_out(i))

       ! The following commands may fail. That is OK. It should just
       ! mean that the attribute is not defined in the input file.

       ncerr = nf90_copy_att(ncid_in, varid_in(i), "long_name",&
            & ncid_out, varid_out(i))
       call handle_err_copy_att("long_name")

       ncerr = nf90_copy_att(ncid_in, varid_in(i), "units", ncid_out,&
            & varid_out(i))
       call handle_err_copy_att("units")

       ncerr = nf90_copy_att(ncid_in, varid_in(i), "standard_name", ncid_out,&
            & varid_out(i))
       call handle_err_copy_att("standard_name")
    end do

    ! Global attributes:
    call nf95_copy_att(ncid_in, nf90_global, "Conventions", ncid_out, &
         nf90_global)
    call nf95_copy_att(ncid_in, nf90_global, "title", ncid_out, nf90_global)
    call nf95_copy_att(ncid_in, nf90_global, "source", ncid_out, nf90_global)
    call nf95_put_att(ncid_out, nf90_global, "comment", "Regridded for LMDZ")

    call nf95_enddef(ncid_out)

    ! Write one of the coordinate variables:
    call nf95_put_var(ncid_out, varid_rlatu, rlatu(jjm+1:1:-1) / pi * 180.)
    ! (convert from rad to degrees and sort in increasing order)

  contains

    subroutine handle_err_copy_att(att_name)

      use netcdf, only: nf90_noerr, nf90_strerror

      character(len=*), intent(in):: att_name

      !----------------------------------------

      if (ncerr /= nf90_noerr) then
         print *, "prepare_out " // trim(name_out(i)) &
              // " nf90_copy_att " // att_name // " -- " &
              // trim(nf90_strerror(ncerr))
      end if

    end subroutine handle_err_copy_att

  end subroutine prepare_out

end module regr_lat_time_coefoz_m
