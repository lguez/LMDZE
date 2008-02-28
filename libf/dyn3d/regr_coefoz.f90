module regr_coefoz_m

  ! This module is clean: no C preprocessor directive, no include line.
  ! Author: Lionel GUEZ

  implicit none

  private
  public regr_coefoz

contains

  subroutine regr_coefoz

    ! "regr_coefoz" stands for "regrid coefficients ozone".

    ! This procedure reads from a file parameters for ozone
    ! chemistry and regrids them for LMDZ.
    ! The input fields depend on time, pressure level and
    ! latitude.
    ! We assume that the input fields are step functions
    ! of latitude.
    ! Horizontal regridding is made by averaging on latitude, with a
    ! cosine of latitude factor.
    ! The target horizontal LMDZ grid is the "scalar" grid: "rlonv", "rlatu".
    ! The values of "rlatu" are taken to be the centers of intervals.
    ! The target vertical LMDZ grid is the grid of mid-layers.
    ! The input latitude values are different from the latitude values
    ! of the LMDZ "scalar" grid.
    ! The input data does not depend on longitude, but the pressure
    ! at LMDZ mid-layers does.
    ! Therefore, the values on the LMDZ grid do depend on longitude.
    ! The regridded fields are written to a file.

    ! We assume that in the input file the latitude is in degrees
    ! and the pressure level is in hPa, and that both are strictly
    ! monotonic (as all NetCDF coordinate variables should be).

    use dimens_m, only: iim, jjm, llm
    use comgeom, only: rlatv
    use comconst, only: pi
    use pressure_m, only: p3d, pls
    use regr1_step_av_m, only: regr1_step_av
    use regr1_lint_m, only: regr1_lint
    use netcdf95, only: nf95_open, nf90_nowrite, coordin, nf95_close, &
         nf95_inq_varid, handle_err, nf90_put_var, nf90_get_var

    ! Variables local to the procedure:

    integer ncid_in, ncid_out ! NetCDF IDs for input and output files

    integer n_plev ! number of pressure levels in the input data
    integer n_lat! number of latitudes in the input data

    real, pointer:: latitude(:)
    ! (of input data, converted to rad, sorted in strictly increasing order)

    real, allocatable:: lat_in_edg(:)
    ! (edges of latitude intervals for input data, in rad, in strictly
    ! increasing order)

    real, pointer:: plev(:)
    ! (pressure level of input data, converted to Pa, sorted
    ! in strictly increasing order)

    real, allocatable:: press_in_edg(:)
    ! (edges of pressure intervals for input data, in Pa, in strictly
    ! increasing order)

    logical decr_lat ! decreasing latitude in the input file
    logical decr_plev ! decreasing pressure level in the input file

    real, allocatable:: o3_par_in(:, :, :)
    ! (ozone parameter from the input file)
    ! ("o3_par_in(j, l, month)" is at latitude "latitude(j)" and pressure
    ! level "plev(l)". "month" is between 1 and 12.)

    real o3_par_out(iim + 1, jjm + 1, llm, 12)
    ! (ozone parameter adapted to the LMDZ grid)
    ! (Last dimension is month number.
    ! "o3_par_out(i, j, l, month)" is at longitude "rlonv(i)", latitude
    ! "rlatu(j)" and pressure level "pls(i, j, l)")

    integer i, j
    integer i_v ! index of ozone parameter
    integer, parameter:: n_o3_param = 8 ! number of ozone parameters

    character(len=11) name_in(n_o3_param)
    ! (name of NetCDF variable in the input file)

    character(len=9) name_out(n_o3_param)
    ! (name of NetCDF variable in the output file)

    logical:: stepav_choice(n_o3_param) = .true.
    ! (vertical regridding by step average, otherwise linear interpolation)

    real top_value(n_o3_param) 
    ! (value at 0 pressure, only used for linear interpolation)

    integer varid_in(n_o3_param), varid_out(n_o3_param), ncerr ! for NetCDF

    real, allocatable:: v_regr_lat(:, :, :) ! (jjm + 1, 0:n_plev, 12)
    ! (mean of a variable "v" over a latitude interval)
    ! First dimension is latitude interval.
    ! The latitude interval for "v_regr_lat(j,:, :)" contains "rlatu(j)".
    ! If "j" is between 2 and "jjm" then the interval is:
    ! [rlatv(j), rlatv(j-1)]
    ! If "j" is 1 or "jjm + 1" then the interval is:
    ! [rlatv(1), pi / 2]
    ! or:
    ! [- pi / 2, rlatv(jjm)]
    ! respectively.
    ! "v_regr_lat(:, l, :)" is for pressure interval
    ! "[press_in_edg(l), press_in_edg(l+1)]" or pressure level "plev(l)",
    ! depending on the type of vertical regridding, step average or linear
    ! interpolation.
    ! Last dimension is month number.)

    !---------------------------------

    print *, "Call sequence information: regr_coefoz"

    ! Details for each ozone parameter:
    i_v = 0

    i_v = i_v + 1
    name_in(i_v) = "P_net"
    name_out(i_v) = "P_net_Mob"

    i_v = i_v + 1
    name_in(i_v) = "a2"
    name_out(i_v) = "a2"

    i_v = i_v + 1
    name_in(i_v) = "r"
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
    stepav_choice(i_v) = .false.
    top_value(i_v) = 0.

    i_v = i_v + 1
    name_in(i_v) = "R_Het"
    name_out(i_v) = "R_Het"

    call nf95_open("coefoz_v2_3.nc", nf90_nowrite, ncid_in)

    ! Get coordinates from the input file:

    latitude => coordin(ncid_in, "latitude")
    ! Convert from degrees to rad, because "rlatv" is in rad:
    latitude = latitude / 180. * pi
    n_lat = size(latitude)
    ! We need to supply the latitudes to "stepav" in
    ! increasing order, so invert order if necessary:
    decr_lat = latitude(1) > latitude(n_lat)
    if (decr_lat) latitude = latitude(n_lat:1:-1)

    ! Compute edges of latitude intervals:
    allocate(lat_in_edg(n_lat + 1))
    lat_in_edg(1) = - pi / 2
    forall (j = 2:n_lat) lat_in_edg(j) = (latitude(j - 1) + latitude(j)) / 2
    lat_in_edg(n_lat + 1) = pi / 2
    deallocate(latitude) ! pointer

    plev => coordin(ncid_in, "plev")
    ! Convert from hPa to Pa because "p3d" and "pls" are in Pa:
    plev = plev * 100.
    n_plev = size(plev)
    ! We need to supply the pressure levels to "stepav" in
    ! increasing order, so invert order if necessary:
    decr_plev = plev(1) > plev(n_plev)
    if (decr_plev) plev = plev(n_plev:1:-1)

    ! Compute edges of pressure intervals:
    allocate(press_in_edg(n_plev + 1))
    press_in_edg(1) = 0.
    ! We choose edges halfway in logarithm:
    forall (j = 2:n_plev) press_in_edg(j) = sqrt(plev(j - 1) * plev(j))
    press_in_edg(n_plev + 1) = huge(0.)
    ! (infinity, but any value guaranteed to be greater than the
    ! surface pressure would do)

    ! Get the IDs of ozone parameters in the input file:
    do i_v = 1, n_o3_param
       call nf95_inq_varid(ncid_in, trim(name_in(i_v)), varid_in(i_v))
    end do

    call prepare_out(ncid_in, varid_in, name_out, ncid_out, varid_out)
    allocate(o3_par_in(n_lat, n_plev, 12), v_regr_lat(jjm + 1, 0:n_plev, 12))

    do i_v = 1, n_o3_param
       ! Process ozone parameter "name_in(i_v)"

       ncerr = nf90_get_var(ncid_in, varid_in(i_v), o3_par_in)
       call handle_err("nf90_get_var", ncerr, ncid_in)

       if (decr_lat) o3_par_in = o3_par_in(n_lat:1:-1, :, :)
       if (decr_plev) o3_par_in = o3_par_in(:, n_plev:1:-1, :)

       ! Regrid in latitude:
       ! We average with respect to sine of latitude, which is
       ! equivalent to weighting by cosine of latitude:
       v_regr_lat(jjm+1:1:-1, 1:, :) = regr1_step_av(o3_par_in, &
            xs=sin(lat_in_edg), xt=sin((/- pi / 2, rlatv(jjm:1:-1), pi / 2/)))
       ! (invert order of indices in "v_regr_lat" because "rlatu" is
       ! decreasing)

       ! Regrid in pressure at each horizontal position:

       if (stepav_choice(i_v)) then
          ! Regrid in pressure by averaging a step function of pressure

          ! Poles ("p3d" does not depend on longitude):
          do j = 1, jjm + 1, jjm
             ! Average on pressure, only for first longitude:
             o3_par_out(1, j, llm:1:-1, :) &
                  = regr1_step_av(v_regr_lat(j, 1:, :), press_in_edg, &
                  p3d(1, j, llm+1:1:-1))
             ! (invert order of indices because "p3d" is decreasing)
          end do

          ! Latitudes other than poles ("p3d" depends on longitude):
          do j = 2, jjm
             ! Average on pressure at each longitude:
             do i = 1, iim
                o3_par_out(i, j, llm:1:-1, :) &
                     = regr1_step_av(v_regr_lat(j, 1:, :), press_in_edg, &
                     p3d(i, j, llm+1:1:-1))
                ! (invert order of indices because "p3d" is decreasing)
             end do
          end do
       else
          ! Regrid in pressure by linear interpolation

          ! Complete "v_regr_lat" with the value at 0 pressure:
          v_regr_lat(:, 0, :) = top_value(i_v)

          ! Poles ("pls" does not depend on longitude):
          do j = 1, jjm + 1, jjm
             ! Interpolate in pressure, only for first longitude:
             o3_par_out(1, j, llm:1:-1, :) = regr1_lint(v_regr_lat(j, :, :), &
                  xs=(/0., plev/), xt=pls(1, j, llm:1:-1))
             ! (invert order of indices because "pls" is decreasing)
          end do

          ! Latitudes other than poles ("pls" depends on longitude):
          do j = 2, jjm
             ! Average on pressure at each longitude:
             do i = 1, iim
                o3_par_out(i, j, llm:1:-1, :) &
                     = regr1_lint(v_regr_lat(j, :, :), xs=(/0., plev/), &
                     xt=pls(i, j, llm:1:-1))
                ! (invert order of indices because "pls" is decreasing)
             end do
          end do
       end if

       ! Duplicate pole values on all longitudes:
       o3_par_out(2:, 1, :, :) &
            = spread(o3_par_out(1, 1, :, :), dim=1, ncopies=iim)
       o3_par_out(2:, jjm + 1, :, :) &
            = spread(o3_par_out(1, jjm + 1, :, :), dim=1, ncopies=iim)

       ! Duplicate first longitude to last longitude:
       o3_par_out(iim + 1, 2:jjm, :, :) = o3_par_out(1, 2:jjm, :, :)

       ! Write to file:

       ncerr = nf90_put_var(ncid_out, varid_out(i_v), &
            o3_par_out(:,jjm+1:1:-1, :, :))
       ! (The order of "rlatu" is inverted in the output file)
       call handle_err("nf90_put_var", ncerr, ncid_out) 
    end do

    deallocate(plev) ! pointer
    call nf95_close(ncid_out)
    call nf95_close(ncid_in)

  end subroutine regr_coefoz

  !********************************************

  subroutine prepare_out(ncid_in, varid_in, name_out, ncid_out, varid_out)

    ! This subroutine creates the NetCDF output file, defines
    ! dimensions and variables, and writes coordinate variables and "pls".

    use dimens_m, only: iim, jjm, llm
    use comgeom, only: rlatu, rlonv
    use pressure_m, only: pls
    use comconst, only: pi
    use nrutil, only: assert_eq

    use netcdf95, only: nf95_create, nf90_clobber, nf95_def_dim, &
         nf95_def_var, nf90_float, nf90_int, nf95_put_att, nf95_enddef, &
         nf90_put_var, handle_err, nf90_copy_att, nf95_copy_att, nf90_global

    integer, intent(in):: ncid_in, varid_in(:)
    character(len=*), intent(in):: name_out(:) ! of NetCDF variables
    integer, intent(out):: ncid_out, varid_out(:)

    ! Variables local to the procedure:
    integer ncerr
    integer dimid_rlonv, dimid_rlatu, dimid_sigs, dimid_month
    integer varid_rlonv, varid_rlatu, varid_sigs, varid_month
    integer varid_pls
    integer i, n_o3_param

    !---------------------------

    print *, "Call sequence information: prepare_out"
    n_o3_param = assert_eq(size(varid_in), size(varid_out), &
         size(name_out), "prepare_out") 

    call nf95_create("coefoz_LMDZ.nc", nf90_clobber, ncid_out)

    ! Dimensions:
    call nf95_def_dim(ncid_out, "month", 12, dimid_month)
    call nf95_def_dim(ncid_out, "sigs", llm, dimid_sigs)
    call nf95_def_dim(ncid_out, "rlatu", jjm + 1, dimid_rlatu)
    call nf95_def_dim(ncid_out, "rlonv", iim + 1, dimid_rlonv)

    ! Coordinate variables:

    call nf95_def_var(ncid_out, "month", nf90_float, dimid_month, varid_month)
    call nf95_put_att(ncid_out, varid_month, "units", &
         "calendar_month as %m.%f")
    call nf95_put_att(ncid_out, varid_month, "long_name", "seasonal phase")

    call nf95_def_var(ncid_out, "sigs", nf90_int, dimid_sigs, varid_sigs)
    call nf95_put_att(ncid_out, varid_sigs, "long_name", "s-layer index")

    call nf95_def_var(ncid_out, "rlatu", nf90_float, dimid_rlatu, varid_rlatu)
    call nf95_put_att(ncid_out, varid_rlatu, "units", "degrees_north")
    call nf95_put_att(ncid_out, varid_rlatu, "long_name", "latitude")

    call nf95_def_var(ncid_out, "rlonv", nf90_float, dimid_rlonv, varid_rlonv)
    call nf95_put_att(ncid_out, varid_rlonv, "units", "degrees_east")
    call nf95_put_att(ncid_out, varid_rlonv, "long_name", "longitude")

    ! Primary variables:

    call nf95_def_var(ncid_out, "pls", nf90_float, &
         (/dimid_rlonv, dimid_rlatu, dimid_sigs/), varid_pls)
    call nf95_put_att(ncid_out, varid_pls, "units", "millibar")
    call nf95_put_att(ncid_out, varid_pls, "long_name", &
         "pressure at LMDZ mid-layers")

    do i = 1, n_o3_param
       call nf95_def_var(ncid_out, name_out(i), nf90_float, &
            (/dimid_rlonv, dimid_rlatu, dimid_sigs, dimid_month/),&
            & varid_out(i))
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
    call nf95_copy_att(ncid_in, nf90_global, "title", ncid_out, nf90_global)
    call nf95_copy_att(ncid_in, nf90_global, "source", ncid_out, nf90_global)
    call nf95_put_att(ncid_out, nf90_global, "comment", "Regridded for LMDZ")

    call nf95_enddef(ncid_out)

    ! Coordinate variables:

    ncerr = nf90_put_var(ncid_out, varid_month, (/(i + 0.5, i = 1, 12)/))
    call handle_err("nf90_put_var", ncerr, ncid_out)

    ncerr = nf90_put_var(ncid_out, varid_sigs, (/(i, i = 1, llm)/))
    call handle_err("nf90_put_var", ncerr, ncid_out)

    ncerr = nf90_put_var(ncid_out, varid_rlatu, rlatu(jjm+1:1:-1) / pi * 180.)
    ! (convert from rad to degrees and sort in increasing order)
    call handle_err("nf90_put_var", ncerr, ncid_out)

    ncerr = nf90_put_var(ncid_out, varid_rlonv, rlonv / pi * 180.)
    ! (convert from rad to degrees)
    call handle_err("nf90_put_var", ncerr, ncid_out)

    ! Primary variable:

    ncerr = nf90_put_var(ncid_out, varid_pls, pls(:, jjm+1:1:-1, :) / 100.)
    ! (convert from Pa to hPa)
    call handle_err("nf90_put_var", ncerr, ncid_out)

  contains

    subroutine handle_err_copy_att(att_name)

      use netcdf95, only: nf90_noerr, nf90_strerror

      character(len=*), intent(in):: att_name

      !----------------------------------------

      if (ncerr /= nf90_noerr) then
         print *, "prepare_out " // trim(name_out(i)) &
              // " nf90_copy_att " // att_name // " -- " &
              // trim(nf90_strerror(ncerr))
      end if

    end subroutine handle_err_copy_att

  end subroutine prepare_out

end module regr_coefoz_m
