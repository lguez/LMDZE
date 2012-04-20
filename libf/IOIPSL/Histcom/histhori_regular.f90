module histhori_regular_m

  implicit none

contains

  SUBROUTINE histhori_regular(pfileid, pim, plon, pjm, plat, phname, phtitle, &
       phid)

    ! This subroutine is made to declare a new horizontale grid.
    ! It has to have the same number of points as
    ! the original Thus in this we routine we will only
    ! add two variable (longitude and latitude).
    ! Any variable in the file can thus point to this pair
    ! through an attribute. This routine is very usefull
    ! to allow staggered grids.

    ! INPUT

    ! pfileid: The id of the file to which the grid should be added
    ! pim: Size in the longitude direction
    ! plon: The longitudes
    ! pjm: Size in the latitude direction
    ! plat: The latitudes
    ! phname: The name of grid
    ! phtitle: The title of the grid

    ! OUTPUT

    ! phid: Id of the created grid

    ! We assume that the grid is rectilinear.

    USE errioipsl, ONLY: histerr
    USE histcom_var, ONLY: full_size, hax_name, nb_hax, ncdf_ids, &
         slab_ori, slab_sz, xid, yid
    USE netcdf, ONLY: nf90_def_var, nf90_enddef, nf90_float, &
         nf90_put_att, nf90_put_var, nf90_redef

    INTEGER, INTENT (IN):: pfileid, pim, pjm
    REAL, INTENT (IN), DIMENSION (pim, pjm):: plon, plat
    CHARACTER (len=*), INTENT (IN):: phname, phtitle
    INTEGER, INTENT (OUT):: phid

    CHARACTER (len=25):: lon_name, lat_name
    CHARACTER (len=80):: tmp_title, tmp_name
    INTEGER:: ndim
    INTEGER, DIMENSION (2):: dims(2)
    INTEGER:: nlonid, nlatid
    INTEGER:: orix, oriy, par_szx, par_szy
    INTEGER:: iret, ncid

    !---------------------------------------------------------------------

    ! 1.0 Check that all fits in the buffers

    IF ((pim/=full_size(pfileid, 1)) .OR. (pjm/=full_size(pfileid, 2))) THEN
       CALL histerr(3, 'histhori', &
            'The new horizontal grid does not have the same size', &
            'as the one provided to histbeg. This is not yet ', &
            'possible in the hist package.')
    END IF

    ! 1.1 Create all the variables needed

    ncid = ncdf_ids(pfileid)

    ndim = 2
    dims(1:2) = (/ xid(pfileid), yid(pfileid) /)

    tmp_name = phname
    IF (nb_hax(pfileid)==0) THEN
       lon_name = 'lon'
       lat_name = 'lat'
    ELSE
       lon_name = 'lon_' // trim(tmp_name)
       lat_name = 'lat_' // trim(tmp_name)
    END IF

    ! 1.2 Save the informations

    phid = nb_hax(pfileid) + 1
    nb_hax(pfileid) = phid

    hax_name(pfileid, phid, 1:2) = (/ lon_name, lat_name/)
    tmp_title = phtitle

    ! 2.0 Longitude

    ndim = 1
    dims(1:1) = (/ xid(pfileid) /)

    iret = nf90_def_var(ncid, lon_name, nf90_float, dims(1:ndim), nlonid)
    iret = nf90_put_att(ncid, nlonid, 'units', 'degrees_east')
    iret = nf90_put_att(ncid, nlonid, 'valid_min', real(minval(plon)))
    iret = nf90_put_att(ncid, nlonid, 'valid_max', real(maxval(plon)))
    iret = nf90_put_att(ncid, nlonid, 'long_name', 'Longitude')
    iret = nf90_put_att(ncid, nlonid, 'nav_model', trim(tmp_title))

    ! 3.0 Latitude

    ndim = 1
    dims(1:1) = (/ yid(pfileid) /)

    iret = nf90_def_var(ncid, lat_name, nf90_float, dims(1:ndim), nlatid)
    iret = nf90_put_att(ncid, nlatid, 'units', 'degrees_north')
    iret = nf90_put_att(ncid, nlatid, 'valid_min', real(minval(plat)))
    iret = nf90_put_att(ncid, nlatid, 'valid_max', real(maxval(plat)))
    iret = nf90_put_att(ncid, nlatid, 'long_name', 'Latitude')
    iret = nf90_put_att(ncid, nlatid, 'nav_model', trim(tmp_title))

    iret = nf90_enddef(ncid)

    ! 4.0 storing the geographical coordinates

    orix = slab_ori(pfileid, 1)
    oriy = slab_ori(pfileid, 2)
    par_szx = slab_sz(pfileid, 1)
    par_szy = slab_sz(pfileid, 2)

    ! Transfer the longitude

    iret = nf90_put_var(ncid, nlonid, plon(1:par_szx, 1))

    ! Transfer the latitude

    iret = nf90_put_var(ncid, nlatid, plat(1, 1:par_szy))

    iret = nf90_redef(ncid)

  END SUBROUTINE histhori_regular

end module histhori_regular_m
