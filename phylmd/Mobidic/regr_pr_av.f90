module regr_pr_av_m

  implicit none

contains

  subroutine regr_pr_av(ncid, name, julien, paprs, v3)

    ! "regr_pr_av" stands for "regrid pressure averaging".

    ! This procedure reads a 2D latitude-pressure field from a NetCDF
    ! file, at a given day, regrids this field in pressure to the LMDZ
    ! vertical grid and packs it to the LMDZ horizontal "physics"
    ! grid. 

    ! We assume that, in the input file, the field has 3 dimensions:
    ! latitude, pressure, julian day.

    ! We assume that the input field is already on the LMDZ "rlatu"
    ! latitudes, except that latitudes are in ascending order in the
    ! input file.

    ! The target vertical LMDZ grid is the grid of layer boundaries.
    ! Regridding in pressure is done by averaging a step function of pressure.

    use dimensions, only: iim, jjm, llm
    use dimphy, only: klon
    use grid_change, only: gr_dyn_phy
    use netcdf95, only: nf95_inq_varid, nf95_get_var
    use jumble, only: assert
    use numer_rec_95, only: regr1_step_av
    use press_coefoz_m, only: press_in_edg

    integer, intent(in):: ncid ! NetCDF ID of the file
    character(len=*), intent(in):: name ! of the NetCDF variable
    integer, intent(in):: julien ! jour julien, 1 <= julien <= 360

    real, intent(in):: paprs(:, :) ! (klon, llm + 1)
    ! (pression pour chaque inter-couche, en Pa)

    real, intent(out):: v3(:, :) ! (klon, llm)
    ! regridded field on the partial "physics" grid
    ! "v3(i, k)" is at longitude "xlon(i)", latitude "xlat(i)", in
    ! layer "k".

    ! Variables local to the procedure:

    integer varid ! for NetCDF

    real v1(jjm + 1, size(press_in_edg) - 1)
    ! input field at day "julien"
    ! "v1(j, k)" is at latitude "rlatu(j)" and for
    ! pressure interval "[press_in_edg(k), press_in_edg(k+1)]".

    real v2(klon, size(press_in_edg) - 1)
    ! Field on the "physics" horizontal grid. "v2(i, k)" is at
    ! longitude "xlon(i)", latitude "xlat(i)" and for pressure
    ! interval "[press_in_edg(k), press_in_edg(k+1)]".)

    integer i

    !--------------------------------------------

    call assert(shape(v3) == (/klon, llm/), "regr_pr_av klon llm")
    call assert(shape(paprs) == (/klon, llm+1/), "regr_pr_av paprs")

    call nf95_inq_varid(ncid, name, varid)

    ! Get data at the right day from the input file:
    call nf95_get_var(ncid, varid, v1, start=(/1, 1, julien/))
    ! Latitudes are in ascending order in the input file while
    ! "rlatu" is in descending order so we need to invert order:
    v1 = v1(jjm+1:1:-1, :)

    v2 = gr_dyn_phy(spread(v1, dim = 1, ncopies = iim + 1))

    ! Regrid in pressure at each horizontal position:
    do i = 1, klon
       v3(i, llm:1:-1) = regr1_step_av(v2(i, :), press_in_edg, &
            paprs(i, llm+1:1:-1))
       ! (invert order of indices because "paprs" is in descending order)
    end do

  end subroutine regr_pr_av

end module regr_pr_av_m
