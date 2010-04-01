MODULE histcom

  ! From histcom.f90, version 2.1 2004/04/21 09:27:10

  ! Some confusing vocabulary in this code !

  ! A REGULAR grid is a grid which is i, j indices
  ! and thus it is stored in a 2D matrix.
  ! This is opposed to an IRREGULAR grid which is only in a vector
  ! and where we do not know which neighbors we have.
  ! As a consequence we need the bounds for each grid-cell.

  ! A RECTILINEAR grid is a special case of a regular grid
  ! in which all longitudes for i constant are equal
  ! and all latitudes for j constant.
  ! In other words we do not need the full 2D matrix
  ! to describe the grid, just two vectors.

  IMPLICIT NONE

  PRIVATE
  PUBLIC:: histbeg_totreg, histdef, histhori_regular, histvert, histend, &
       histclo, histsync

CONTAINS

  SUBROUTINE histbeg_totreg(pfilename, plon_1d, plat_1d, par_orix, par_szx, &
       par_oriy, par_szy, pitau0, pdate0, pdeltat, phoriid, pfileid)

    ! The user provides "plon_1d" and "plat_1d" as vectors. Obviously
    ! this can only be used for very regular grids.
    ! This subroutine initializes a netcdf file and returns the ID.
    ! It will set up the geographical space on which the data will be
    ! stored and offers the possibility of seting a zoom.
    ! It also gets the global parameters into the I/O subsystem.

    ! INPUT

    ! pfilename: Name of the netcdf file to be created
    ! pim: Size of arrays in longitude direction
    ! plon_1d: Coordinates of points in longitude
    ! pjm: Size of arrays in latitude direction
    ! plat_1d: Coordinates of points in latitude

    ! The next 4 arguments allow to define a horizontal zoom
    ! for this file. It is assumed that all variables to come
    ! have the same index space. This can not be assumed for
    ! the z axis and thus we define the zoom in histdef.

    ! par_orix: Origin of the slab of data within the X axis (pim)
    ! par_szx: Size of the slab of data in X
    ! par_oriy: Origin of the slab of data within the Y axis (pjm)
    ! par_szy: Size of the slab of data in Y

    ! pitau0: time step at which the history tape starts
    ! pdate0: The Julian date at which the itau was equal to 0
    ! pdeltat: Time step in seconds. Time step of the counter itau
    !             used in histwrt for instance

    ! OUTPUT

    ! phoriid: ID of the horizontal grid
    ! pfileid: ID of the netcdf file

    ! We assume the grid is rectilinear.

    USE ioipslmpp, ONLY: ioipslmpp_file
    USE errioipsl, ONLY: histerr
    USE histcom_var, ONLY: assc_file, date0, deltat, full_size, itau0, &
         lock_modname, model_name, nb_files, nb_files_max, nb_hax, nb_tax, &
         nb_var, nb_zax, ncdf_ids, regular, slab_ori, slab_sz, xid, yid, zoom
    USE netcdf, ONLY: nf90_clobber, nf90_create, nf90_def_dim, &
         nf90_global, nf90_put_att

    CHARACTER (len=*), INTENT (IN):: pfilename
    REAL, DIMENSION (:), INTENT (IN):: plon_1d
    REAL, DIMENSION (:), INTENT (IN):: plat_1d
    INTEGER, INTENT (IN):: par_orix, par_szx, par_oriy, par_szy
    INTEGER, INTENT (IN):: pitau0
    REAL, INTENT (IN):: pdate0, pdeltat
    INTEGER, INTENT (OUT):: pfileid, phoriid

    ! Variables local to the procedure:
    REAL, DIMENSION (size(plon_1d), size(plat_1d)):: plon, plat
    INTEGER:: pim, pjm
    INTEGER:: ncid, iret
    INTEGER:: lengf, lenga
    CHARACTER (len=120):: file, tfile

    !---------------------------------------------------------------------

    pim = size(plon_1d)
    pjm = size(plat_1d)

    plon = spread(plon_1d, 2, pjm)
    plat = spread(plat_1d, 1, pim)

    nb_files = nb_files + 1
    pfileid = nb_files

    ! 1.0 Transfering into the common for future use

    itau0(pfileid) = pitau0
    date0(pfileid) = pdate0
    deltat(pfileid) = pdeltat

    ! 2.0 Initializes all variables for this file

    IF (nb_files>nb_files_max) THEN
       CALL histerr(3, 'histbeg', &
            'Table of files too small. You should increase nb_files_max', &
            'in M_HISTCOM.f90 in order to accomodate all these files', ' ')
    END IF

    nb_var(pfileid) = 0
    nb_tax(pfileid) = 0
    nb_hax(pfileid) = 0
    nb_zax(pfileid) = 0

    slab_ori(pfileid, 1:2) = (/ par_orix, par_oriy/)
    slab_sz(pfileid, 1:2) = (/ par_szx, par_szy/)

    ! 3.0 Opening netcdf file and defining dimensions

    tfile = pfilename
    lengf = len_trim(tfile)
    IF (tfile(lengf-2:lengf)/='.nc') THEN
       file = tfile(1:lengf) // '.nc'
    ELSE
       file = tfile(1:lengf)
    END IF

    ! Add PE number in file name on MPP

    CALL ioipslmpp_file(file)

    ! Keep track of the name of the files opened

    lengf = len_trim(file)
    lenga = len_trim(assc_file)
    IF (nb_files==1) THEN
       assc_file = file(1:lengf)
    ELSE IF ((lenga+lengf)<500) THEN
       assc_file = assc_file(1:lenga) // ' ' // file(1:lengf)
    ELSE IF (((lenga+7)<500) .AND. (index(assc_file(1:lenga), &
         'et.al.')<1)) THEN
       assc_file = assc_file(1:lenga) // ' et.al.'
    ELSE
       CALL histerr(2, 'histbeg', &
            'The file names do not fit into the associate_file attribute.', &
            'Use shorter names if you wish to keep the information.', ' ')
    END IF

    iret = nf90_create(file, nf90_clobber, ncid)
    iret = nf90_def_dim(ncid, 'lon', par_szx, xid(nb_files))
    iret = nf90_def_dim(ncid, 'lat', par_szy, yid(nb_files))

    ! 4.0 Declaring the geographical coordinates and other attributes

    ! 4.3 Global attributes

    iret = nf90_put_att(ncid, nf90_global, 'Conventions', 'GDT 1.3')
    iret = nf90_put_att(ncid, nf90_global, 'file_name', trim(file))
    iret = nf90_put_att(ncid, nf90_global, 'production', trim(model_name))
    lock_modname = .TRUE.

    ! 5.0 Saving some important information on this file in the common

    ncdf_ids(pfileid) = ncid
    full_size(pfileid, 1:2) = (/ pim, pjm/)

    ! 6.0 storing the geographical coordinates

    IF ((pim/=par_szx) .OR. (pjm/=par_szy)) zoom(pfileid) = .TRUE.
    regular(pfileid) = .TRUE.

    CALL histhori_regular(pfileid, pim, plon, pjm, plat, ' ', 'Default grid', &
         phoriid)

  END SUBROUTINE histbeg_totreg

  !**********************************

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

  !**********************************

  SUBROUTINE histvert(pfileid, pzaxname, pzaxtitle, pzaxunit, pzsize, pzvalues, &
       pzaxid, pdirect)

    ! This subroutine defines a vertical axis and returns it s id.
    ! It gives the user the possibility to the user to define many
    ! different vertical axes. For each variable defined with histdef a
    ! vertical axis can be specified with by it s ID.

    ! INPUT

    ! pfileid: ID of the file the variable should be archived in
    ! pzaxname: Name of the vertical axis
    ! pzaxtitle: title of the vertical axis
    ! pzaxunit: Units of the vertical axis
    ! pzsize: size of the vertical axis
    ! pzvalues: Coordinate values of the vetical axis

    ! pdirect: is an optional argument which allows to specify the
    !            the positive direction of the axis: up or down.
    ! OUTPUT

    ! pzaxid: Returns the ID of the axis.
    !            Note that this is not the netCDF ID !

    USE stringop, ONLY: find_str, strlowercase
    USE stringop, ONLY: find_str, strlowercase
    USE errioipsl, ONLY: histerr
    USE histcom_var, ONLY: nb_zax, nb_zax_max, ncdf_ids, zax_ids, &
         zax_name, zax_name_length, zax_size
    USE netcdf, ONLY: nf90_def_dim, nf90_def_var, nf90_enddef, &
         nf90_float, nf90_put_att, nf90_put_var, nf90_redef

    INTEGER, INTENT (IN):: pfileid, pzsize
    CHARACTER (len=*), INTENT (IN):: pzaxname, pzaxunit, pzaxtitle
    REAL, INTENT (IN):: pzvalues(pzsize)
    INTEGER, INTENT (OUT):: pzaxid
    CHARACTER (len=*), INTENT (IN), OPTIONAL:: pdirect

    INTEGER:: pos, iv, nb, zdimid, zaxid_tmp
    CHARACTER (len=20):: str20, tab_str20(nb_zax_max)
    INTEGER:: tab_str20_length(nb_zax_max)
    CHARACTER (len=70):: str70, str71, str72
    CHARACTER (len=80):: str80
    CHARACTER (len=20):: direction
    INTEGER:: iret, leng, ncid

    !---------------------------------------------------------------------

    ! 1.0 Verifications:
    !    Do we have enough space for an extra axis ?
    !    Is the name already in use ?

    ! - Direction of axis. Can we get if from the user.
    !   If not we put unknown.

    IF (present(pdirect)) THEN
       direction = trim(pdirect)
       CALL strlowercase(direction)
    ELSE
       direction = 'unknown'
    END IF

    ! Check the consistency of the attribute

    IF ((direction/='unknown') .AND. (direction/='up') .AND. &
         (direction/='down')) THEN
       direction = 'unknown'
       str80 = 'The specified axis was: ' // trim(direction)
       CALL histerr(2, 'histvert', &
            'The specified direction for the vertical axis is not possible.', &
            'it is replaced by: unknown', str80)
    END IF

    IF (nb_zax(pfileid)+1>nb_zax_max) THEN
       CALL histerr(3, 'histvert', &
            'Table of vertical axes too small. You should increase ', &
            'nb_zax_max in M_HISTCOM.f90 in order to accomodate all ', &
            'these variables ')
    END IF

    iv = nb_zax(pfileid)
    IF (iv>1) THEN
       str20 = pzaxname
       nb = iv - 1
       tab_str20(1:nb) = zax_name(pfileid, 1:nb)
       tab_str20_length(1:nb) = zax_name_length(pfileid, 1:nb)
       CALL find_str(nb, tab_str20, tab_str20_length, str20, pos)
    ELSE
       pos = 0
    END IF

    IF (pos>0) THEN
       str70 = 'Vertical axis already exists'
       WRITE (str71, '("Check variable ", a, " in file", I3)') str20, &
            pfileid
       str72 = 'Can also be a wrong file ID in another declaration'
       CALL histerr(3, 'histvert', str70, str71, str72)
    END IF

    iv = nb_zax(pfileid) + 1

    ! 2.0 Add the information to the file

    ncid = ncdf_ids(pfileid)

    leng = min(len_trim(pzaxname), 20)
    iret = nf90_def_dim(ncid, pzaxname(1:leng), pzsize, zaxid_tmp)
    iret = nf90_def_var(ncid, pzaxname(1:leng), nf90_float, zaxid_tmp, zdimid)

    leng = min(len_trim(pzaxunit), 20)
    iret = nf90_put_att(ncid, zdimid, 'units', pzaxunit(1:leng))
    iret = nf90_put_att(ncid, zdimid, 'positive', trim(direction))

    iret = nf90_put_att(ncid, zdimid, 'valid_min', real(minval( &
         pzvalues(1:pzsize))))
    iret = nf90_put_att(ncid, zdimid, 'valid_max', real(maxval( &
         pzvalues(1:pzsize))))

    leng = min(len_trim(pzaxname), 20)
    iret = nf90_put_att(ncid, zdimid, 'title', pzaxname(1:leng))
    leng = min(len_trim(pzaxtitle), 80)
    iret = nf90_put_att(ncid, zdimid, 'long_name', pzaxtitle(1:leng))

    iret = nf90_enddef(ncid)

    iret = nf90_put_var(ncid, zdimid, pzvalues(1:pzsize))

    iret = nf90_redef(ncid)

    ! 3.0 add the information to the common

    nb_zax(pfileid) = iv
    zax_size(pfileid, iv) = pzsize
    zax_name(pfileid, iv) = pzaxname
    zax_name_length(pfileid, iv) = len_trim(pzaxname)
    zax_ids(pfileid, iv) = zaxid_tmp
    pzaxid = iv

  END SUBROUTINE histvert

  !**********************************

  SUBROUTINE histdef(pfileid, pvarname, ptitle, punit, pxsize, pysize, phoriid, &
       pzsize, par_oriz, par_szz, pzid, popp, pfreq_opp, pfreq_wrt)

    ! With this subroutine each variable to be archived on the history
    ! tape should be declared.

    ! It gives the user the choise of operation
    ! to be performed on the variables, the frequency of this operation
    ! and finaly the frequency of the archiving.

    USE stringop, ONLY: find_str
    USE mathelp, ONLY: buildop
    USE errioipsl, ONLY: histerr
    USE histcom_var, ONLY: buff_pos, deltat, freq_opp, freq_wrt, fullop, &
         full_size, itau0, last_opp, last_opp_chk, last_wrt, last_wrt_chk, &
         missing_val, name, name_length, nbopp, nbopp_max, nb_hax, nb_opp, &
         nb_tax, nb_var, nb_var_max, nb_wrt, nb_zax, point, scal, scsize, &
         slab_ori, slab_sz, sopps, tax_last, tax_name, tax_name_length, &
         title, topp, unit_name, var_axid, var_haxid, var_zaxid, zax_name, &
         zax_size, zorig, zsize
    USE calendar, ONLY: ioget_calendar

    INTEGER, INTENT (IN):: pfileid
    ! (ID of the file the variable should be archived in)

    CHARACTER (len=*), INTENT (IN):: pvarname
    ! (Name of the variable, short and easy to remember)

    CHARACTER (len=*), INTENT (IN):: ptitle ! Full name of the variable
    CHARACTER (len=*), INTENT (IN):: punit ! Units of the variable

    ! The next 3 arguments give the size of that data
    ! that will be passed to histwrite. The zoom will be
    ! done there with the horizontal information obtained
    ! in "histbeg" and the vertical information to follow.
    INTEGER, INTENT (IN):: pxsize, pysize ! Sizes in X and Y directions
    INTEGER, INTENT (IN):: phoriid ! ID of the horizontal axis

    ! The next two arguments give the vertical zoom to use.

    INTEGER, INTENT (IN):: pzsize
    ! (Size in Z direction (If 1 then no axis is declared for this
    ! variable and pzid is not used)

    INTEGER, INTENT (IN):: par_oriz ! Off set of the zoom
    INTEGER, INTENT (IN):: par_szz ! Size of the zoom

    INTEGER, INTENT (IN):: pzid
    ! (ID of the vertical axis to use. It has to have the size of the zoom.)

    CHARACTER (len=*), INTENT (IN):: popp
    ! Operation to be performed. The following options exist today:
    ! inst: keeps instantaneous values for writting
    ! ave: Computes the average from call between writes

    REAL, INTENT (IN):: pfreq_opp ! Frequency of this operation (in seconds)

    REAL, INTENT (IN):: pfreq_wrt
    ! (Frequency at which the variable should be written, in seconds)

    INTEGER:: iv, i, nb
    CHARACTER (len=70):: str70, str71, str72
    CHARACTER (len=20):: tmp_name
    CHARACTER (len=20):: str20, tab_str20(nb_var_max)
    INTEGER:: tab_str20_length(nb_var_max)
    CHARACTER (len=40):: str40, tab_str40(nb_var_max)
    INTEGER:: tab_str40_length(nb_var_max)
    CHARACTER (len=10):: str10
    CHARACTER (len=80):: tmp_str80
    CHARACTER (len=7):: tmp_topp, tmp_sopp(nbopp_max)
    CHARACTER (len=120):: ex_topps
    REAL:: tmp_scal(nbopp_max), un_an, un_jour, test_fopp, test_fwrt
    INTEGER:: pos, buff_sz

    !---------------------------------------------------------------------
    ex_topps = 'ave, inst, t_min, t_max, t_sum, once, never, l_max, l_min'

    nb_var(pfileid) = nb_var(pfileid) + 1
    iv = nb_var(pfileid)

    IF (iv>nb_var_max) THEN
       CALL histerr(3, 'histdef', &
            'Table of variables too small. You should increase nb_var_max', &
            'in M_HISTCOM.f90 in order to accomodate all these variables', ' ')
    END IF

    ! 1.0 Transfer informations on the variable to the common
    !     and verify that it does not already exist

    IF (iv>1) THEN
       str20 = pvarname
       nb = iv - 1
       tab_str20(1:nb) = name(pfileid, 1:nb)
       tab_str20_length(1:nb) = name_length(pfileid, 1:nb)
       CALL find_str(nb, tab_str20, tab_str20_length, str20, pos)
    ELSE
       pos = 0
    END IF

    IF (pos>0) THEN
       str70 = 'Variable already exists'
       WRITE (str71, '("Check variable  ", a, " in file", I3)') str20, &
            pfileid
       str72 = 'Can also be a wrong file ID in another declaration'
       CALL histerr(3, 'histdef', str70, str71, str72)
    END IF

    name(pfileid, iv) = pvarname
    name_length(pfileid, iv) = len_trim(name(pfileid, iv))
    title(pfileid, iv) = ptitle
    unit_name(pfileid, iv) = punit
    tmp_name = name(pfileid, iv)

    ! 1.1 decode the operations

    fullop(pfileid, iv) = popp
    tmp_str80 = popp
    CALL buildop(tmp_str80, ex_topps, tmp_topp, nbopp_max, missing_val, &
         tmp_sopp, tmp_scal, nbopp(pfileid, iv))

    topp(pfileid, iv) = tmp_topp
    DO i = 1, nbopp(pfileid, iv)
       sopps(pfileid, iv, i) = tmp_sopp(i)
       scal(pfileid, iv, i) = tmp_scal(i)
    END DO

    ! 1.2 If we have an even number of operations
    !     then we need to add identity

    IF (2*int(nbopp(pfileid, iv)/2.0)==nbopp(pfileid, iv)) THEN
       nbopp(pfileid, iv) = nbopp(pfileid, iv) + 1
       sopps(pfileid, iv, nbopp(pfileid, iv)) = 'ident'
       scal(pfileid, iv, nbopp(pfileid, iv)) = missing_val
    END IF

    ! 2.0 Put the size of the variable in the common and check

    scsize(pfileid, iv, :) = (/ pxsize, pysize, pzsize/)

    zorig(pfileid, iv, 1:3) = (/ slab_ori(pfileid, 1), slab_ori(pfileid, 2), &
         par_oriz/)

    zsize(pfileid, iv, 1:3) = (/ slab_sz(pfileid, 1), slab_sz(pfileid, 2), &
         par_szz/)

    ! Is the size of the full array the same as that of the coordinates  ?

    IF ((pxsize>full_size(pfileid, 1)) .OR. (pysize>full_size(pfileid, &
         2))) THEN

       str70 = 'The size of the variable is different ' // &
            'from the one of the coordinates'
       WRITE (str71, '("Size of coordinates:", 2I4)') full_size(pfileid, 1), &
            full_size(pfileid, 2)
       WRITE (str72, '("Size declared for variable ", a, ":", 2I4)') &
            trim(tmp_name), pxsize, pysize
       CALL histerr(3, 'histdef', str70, str71, str72)
    END IF

    ! Is the size of the zoom smaler than the coordinates ?

    IF ((full_size(pfileid, 1)<slab_sz(pfileid, 1)) .OR. (full_size(pfileid, &
         2)<slab_sz(pfileid, 2))) THEN
       str70 = 'Size of variable should be greater or equal &
            &to those of the zoom'
       WRITE (str71, '("Size of XY zoom:", 2I4)') slab_sz(pfileid, 1), &
            slab_sz(pfileid, 1)
       WRITE (str72, '("Size declared for variable ", a, ":", 2I4)') &
            trim(tmp_name), pxsize, pysize
       CALL histerr(3, 'histdef', str70, str71, str72)
    END IF

    ! 2.1 We store the horizontal grid information with minimal
    !     and a fall back onto the default grid

    IF (phoriid>0 .AND. phoriid<=nb_hax(pfileid)) THEN
       var_haxid(pfileid, iv) = phoriid
    ELSE
       var_haxid(pfileid, iv) = 1
       CALL histerr(2, 'histdef', &
            'We use the default grid for variable as an invalide', &
            'ID was provided for variable: ', pvarname)
    END IF

    ! 2.2 Check the vertical coordinates if needed

    IF (par_szz>1) THEN

       ! Does the vertical coordinate exist ?

       IF (pzid>nb_zax(pfileid)) THEN
          WRITE (str70, '("The vertical coordinate chosen for variable ", a)' &
               ) trim(tmp_name)
          str71 = ' Does not exist.'
          CALL histerr(3, 'histdef', str70, str71, ' ')
       END IF

       ! Is the vertical size of the variable equal to that of the axis ?

       IF (par_szz/=zax_size(pfileid, pzid)) THEN
          str20 = zax_name(pfileid, pzid)
          str70 = 'The size of the zoom does not correspond ' // &
               'to the size of the chosen vertical axis'
          WRITE (str71, '("Size of zoom in z:", I4)') par_szz
          WRITE (str72, '("Size declared for axis ", a, ":", I4)') &
               trim(str20), zax_size(pfileid, pzid)
          CALL histerr(3, 'histdef', str70, str71, str72)
       END IF

       ! Is the zoom smaler that the total size of the variable ?

       IF (pzsize<par_szz) THEN
          str20 = zax_name(pfileid, pzid)
          str70 = 'The vertical size of variable ' // &
               'is smaller than that of the zoom.'
          WRITE (str71, '("Declared vertical size of data:", I5)') pzsize
          WRITE (str72, '("Size of zoom for variable ", a, " = ", I5)') &
               trim(tmp_name), par_szz
          CALL histerr(3, 'histdef', str70, str71, str72)
       END IF
       var_zaxid(pfileid, iv) = pzid
    ELSE
       var_zaxid(pfileid, iv) = -99
    END IF

    ! 3.0 Determine the position of the variable in the buffer
    !     If it is instantaneous output then we do not use the buffer

    ! 3.1 We get the size of the arrays histwrite will get and check
    !     that they fit into the tmp_buffer

    buff_sz = zsize(pfileid, iv, 1)*zsize(pfileid, iv, 2)*zsize(pfileid, iv, 3)

    ! 3.2 move the pointer of the buffer array for operation
    !     which need bufferisation

    IF ((trim(tmp_topp)/='inst') .AND. (trim(tmp_topp)/='once') .AND. ( &
         trim(tmp_topp)/='never')) THEN
       point(pfileid, iv) = buff_pos + 1
       buff_pos = buff_pos + buff_sz
    END IF

    ! 4.0 Transfer the frequency of the operations and check
    !     for validity. We have to pay attention to negative values
    !     of the frequency which indicate monthly time-steps.
    !     The strategy is to bring it back to seconds for the tests

    freq_opp(pfileid, iv) = pfreq_opp
    freq_wrt(pfileid, iv) = pfreq_wrt

    CALL ioget_calendar(un_an, un_jour)
    IF (pfreq_opp<0) THEN
       CALL ioget_calendar(un_an)
       test_fopp = pfreq_opp*(-1.)*un_an/12.*un_jour
    ELSE
       test_fopp = pfreq_opp
    END IF
    IF (pfreq_wrt<0) THEN
       CALL ioget_calendar(un_an)
       test_fwrt = pfreq_wrt*(-1.)*un_an/12.*un_jour
    ELSE
       test_fwrt = pfreq_wrt
    END IF

    ! 4.1 Frequency of operations and output should be larger than deltat !

    IF (test_fopp<deltat(pfileid)) THEN
       str70 = 'Frequency of operations should be larger than deltat'
       WRITE (str71, '("It is not the case for variable ", a, ":", F10.4)') &
            trim(tmp_name), pfreq_opp
       str72 = 'PATCH: frequency set to deltat'

       CALL histerr(2, 'histdef', str70, str71, str72)

       freq_opp(pfileid, iv) = deltat(pfileid)
    END IF

    IF (test_fwrt<deltat(pfileid)) THEN
       str70 = 'Frequency of output should be larger than deltat'
       WRITE (str71, '("It is not the case for variable ", a, ":", F10.4)') &
            trim(tmp_name), pfreq_wrt
       str72 = 'PATCH: frequency set to deltat'

       CALL histerr(2, 'histdef', str70, str71, str72)

       freq_wrt(pfileid, iv) = deltat(pfileid)
    END IF

    ! 4.2 First the existence of the operation is tested and then
    !     its compaticility with the choice of frequencies

    IF (trim(tmp_topp)=='inst') THEN
       IF (test_fopp/=test_fwrt) THEN
          str70 = 'For instantaneous output the frequency ' // &
               'of operations and output'
          WRITE (str71, &
               '("should be the same, this was not case for variable ", a)') &
               trim(tmp_name)
          str72 = 'PATCH: The smalest frequency of both is used'
          CALL histerr(2, 'histdef', str70, str71, str72)
          IF (test_fopp<test_fwrt) THEN
             freq_opp(pfileid, iv) = pfreq_opp
             freq_wrt(pfileid, iv) = pfreq_opp
          ELSE
             freq_opp(pfileid, iv) = pfreq_wrt
             freq_wrt(pfileid, iv) = pfreq_wrt
          END IF
       END IF
    ELSE IF (index(ex_topps, trim(tmp_topp))>0) THEN
       IF (test_fopp>test_fwrt) THEN
          str70 = 'For averages the frequency of operations ' // &
               'should be smaller or equal'
          WRITE (str71, &
               '("to that of output. It is not the case for variable ", a)') &
               trim(tmp_name)
          str72 = 'PATCH: The output frequency is used for both'
          CALL histerr(2, 'histdef', str70, str71, str72)
          freq_opp(pfileid, iv) = pfreq_wrt
       END IF
    ELSE
       WRITE (str70, '("Operation on variable ", a, " is unknown")') &
            trim(tmp_name)
       WRITE (str71, '("operation requested is:", a)') tmp_topp
       WRITE (str72, '("File ID:", I3)') pfileid
       CALL histerr(3, 'histdef', str70, str71, str72)
    END IF

    ! 5.0 Initialize other variables of the common

    last_opp(pfileid, iv) = itau0(pfileid)
    ! - freq_opp(pfileid, iv)/2./deltat(pfileid)
    last_wrt(pfileid, iv) = itau0(pfileid)
    ! - freq_wrt(pfileid, iv)/2./deltat(pfileid)
    last_opp_chk(pfileid, iv) = itau0(pfileid)
    ! - freq_opp(pfileid, iv)/2./deltat(pfileid)
    last_wrt_chk(pfileid, iv) = itau0(pfileid)
    ! - freq_wrt(pfileid, iv)/2./deltat(pfileid)
    nb_opp(pfileid, iv) = 0
    nb_wrt(pfileid, iv) = 0

    ! 6.0 Get the time axis for this variable

    IF (freq_wrt(pfileid, iv)>0) THEN
       WRITE (str10, '(I8.8)') int(freq_wrt(pfileid, iv))
       str40 = trim(tmp_topp) // '_' // trim(str10)
    ELSE
       WRITE (str10, '(I2.2, "month")') abs(int(freq_wrt(pfileid, iv)))
       str40 = trim(tmp_topp) // '_' // trim(str10)
    END IF

    DO i = 1, nb_tax(pfileid)
       tab_str40(i) = tax_name(pfileid, i)
       tab_str40_length(i) = tax_name_length(pfileid, i)
    END DO

    CALL find_str(nb_tax(pfileid), tab_str40, tab_str40_length, str40, pos)

    ! No time axis for once, l_max, l_min or never operation

    IF ((trim(tmp_topp)/='once') .AND. (trim(tmp_topp)/='never') .AND. ( &
         trim(tmp_topp)/='l_max') .AND. (trim(tmp_topp)/='l_min')) THEN
       IF (pos<0) THEN
          nb_tax(pfileid) = nb_tax(pfileid) + 1
          tax_name(pfileid, nb_tax(pfileid)) = str40
          tax_name_length(pfileid, nb_tax(pfileid)) = len_trim(str40)
          tax_last(pfileid, nb_tax(pfileid)) = 0
          var_axid(pfileid, iv) = nb_tax(pfileid)
       ELSE
          var_axid(pfileid, iv) = pos
       END IF
    ELSE
       var_axid(pfileid, iv) = -99
    END IF

    ! 7.0 prepare frequence of writing and operation
    !     for never or once operation

    IF ((trim(tmp_topp)=='once') .OR. (trim(tmp_topp)=='never')) THEN
       freq_opp(pfileid, iv) = 0.
       freq_wrt(pfileid, iv) = 0.
    END IF

  END SUBROUTINE histdef

  !**********************************

  SUBROUTINE histend(pfileid)

    ! This subroutine ends the declaration of variables and sets the
    ! time axes in the netcdf file and puts it into write mode.

    ! INPUT

    ! pfileid: ID of the file to be worked on

    USE ioipslmpp, ONLY: ioipslmpp_addatt
    USE errioipsl, ONLY: histerr
    USE histcom_var, ONLY: date0, freq_opp, freq_wrt, fullop, &
         missing_val, name, nb_tax, nb_var, ncdf_ids, ncvar_ids, regular, &
         tax_name, tdimid, tid, title, topp, unit_name, var_axid, var_zaxid, &
         xid, yid, zax_ids, zax_name
    USE calendar, ONLY: ioget_calendar, ju2ymds
    USE netcdf, ONLY: nf90_def_dim, nf90_def_var, nf90_enddef, &
         nf90_float, nf90_put_att, nf90_unlimited

    INTEGER, INTENT (IN):: pfileid

    INTEGER:: ncid, ncvarid
    INTEGER:: iret, ndim, iv, itx, ziv
    INTEGER:: itax
    INTEGER:: dims(4), dim_cnt
    INTEGER:: year, month, day, hours, minutes
    REAL:: sec
    REAL:: rtime0
    CHARACTER (len=20):: tname, tunit
    CHARACTER (len=30):: str30
    CHARACTER (len=80):: ttitle
    CHARACTER (len=120):: assoc
    CHARACTER (len=70):: str70
    CHARACTER (len=3), DIMENSION (12):: cal = (/ 'JAN', 'FEB', 'MAR', &
         'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'/)
    CHARACTER (len=7):: tmp_opp

    !---------------------------------------------------------------------
    ncid = ncdf_ids(pfileid)

    ! 1.0 Create the time axes

    iret = nf90_def_dim(ncid, 'time_counter', nf90_unlimited, tid(pfileid))

    ! 1.1 Define all the time axes needed for this file

    DO itx = 1, nb_tax(pfileid)
       dims(1) = tid(pfileid)
       IF (nb_tax(pfileid)>1) THEN
          str30 = 't_' // tax_name(pfileid, itx)
       ELSE
          str30 = 'time_counter'
       END IF
       iret = nf90_def_var(ncid, str30, nf90_float, dims(1), &
            tdimid(pfileid, itx))

       !   To transform the current itau into a real date and take it
       !   as the origin of the file requires the time counter to change.
       !   Thus it is an operation the user has to ask for.
       !   This function should thus only be re-instated
       !   if there is a ioconf routine to control it.

       ! rtime0 = itau2date(itau0(pfileid), date0(pfileid), deltat(pfileid))
       rtime0 = date0(pfileid)

       CALL ju2ymds(rtime0, year, month, day, sec)

       !   Catch any error induced by a change in calendar !

       IF (year<0) THEN
          year = 2000 + year
       END IF

       hours = int(sec/(60.*60.))
       minutes = int((sec-hours*60.*60.)/60.)
       sec = sec - (hours*60.*60.+minutes*60.)

       WRITE (str70, 7000) year, month, day, hours, minutes, int(sec)
       iret = nf90_put_att(ncid, tdimid(pfileid, itx), 'units', trim(str70))

       CALL ioget_calendar(str30)
       iret = nf90_put_att(ncid, tdimid(pfileid, itx), 'calendar', trim(str30))

       iret = nf90_put_att(ncid, tdimid(pfileid, itx), 'title', 'Time')

       iret = nf90_put_att(ncid, tdimid(pfileid, itx), 'long_name', &
            'Time axis')

       WRITE (str70, 7001) year, cal(month), day, hours, minutes, int(sec)
       iret = nf90_put_att(ncid, tdimid(pfileid, itx), 'time_origin', &
            trim(str70))
    END DO

    ! The formats we need

7000 FORMAT ('seconds since ', I4.4, '-', I2.2, '-', I2.2, ' ', I2.2, ':', I2.2, ':', &
         I2.2)
7001 FORMAT (' ', I4.4, '-', A3, '-', I2.2, ' ', I2.2, ':', I2.2, ':', I2.2)

    ! 2.0 declare the variables

    DO iv = 1, nb_var(pfileid)

       itax = var_axid(pfileid, iv)

       tname = name(pfileid, iv)
       tunit = unit_name(pfileid, iv)
       ttitle = title(pfileid, iv)

       IF (regular(pfileid)) THEN
          dims(1:2) = (/ xid(pfileid), yid(pfileid) /)
          dim_cnt = 2
       ELSE
          dims(1) = xid(pfileid)
          dim_cnt = 1
       END IF

       tmp_opp = topp(pfileid, iv)
       ziv = var_zaxid(pfileid, iv)

       !   2.1 dimension of field

       IF ((trim(tmp_opp)/='never')) THEN
          IF ((trim(tmp_opp)/='once') .AND. (trim( &
               tmp_opp)/='l_max') .AND. (trim(tmp_opp)/='l_min')) THEN
             IF (ziv==-99) THEN
                ndim = dim_cnt + 1
                dims(dim_cnt+1:dim_cnt+2) = (/ tid(pfileid), 0 /)
             ELSE
                ndim = dim_cnt + 2
                dims(dim_cnt+1:dim_cnt+2) = (/ zax_ids(pfileid, ziv), &
                     tid(pfileid) /)
             END IF
          ELSE
             IF (ziv==-99) THEN
                ndim = dim_cnt
                dims(dim_cnt+1:dim_cnt+2) = (/ 0, 0 /)
             ELSE
                ndim = dim_cnt + 1
                dims(dim_cnt+1:dim_cnt+2) = (/ zax_ids(pfileid, ziv), 0 /)
             END IF
          END IF

          iret = nf90_def_var(ncid, trim(tname), nf90_float, dims(1:abs(ndim)), &
               ncvarid)

          ncvar_ids(pfileid, iv) = ncvarid

          iret = nf90_put_att(ncid, ncvarid, 'units', trim(tunit))

          iret = nf90_put_att(ncid, ncvarid, 'missing_value', &
               real(missing_val))
          iret = nf90_put_att(ncid, ncvarid, 'long_name', trim(ttitle))

          iret = nf90_put_att(ncid, ncvarid, 'short_name', trim(tname))

          iret = nf90_put_att(ncid, ncvarid, 'online_operation', trim(fullop( &
               pfileid, iv)))

          SELECT CASE (ndim)
          CASE (-3)
             str30 = 'ZYX'
          CASE (2)
             str30 = 'YX'
          CASE (3)
             str30 = 'TYX'
          CASE (4)
             str30 = 'TZYX'
          CASE DEFAULT
             CALL histerr(3, 'histend', &
                  'less than 2 or more than 4 dimensions are not', &
                  'allowed at this stage', ' ')
          END SELECT

          iret = nf90_put_att(ncid, ncvarid, 'axis', trim(str30))

          assoc = 'nav_lat nav_lon'
          ziv = var_zaxid(pfileid, iv)
          IF (ziv>0) THEN
             str30 = zax_name(pfileid, ziv)
             assoc = trim(str30) // ' ' // trim(assoc)
          END IF

          IF (itax>0) THEN
             IF (nb_tax(pfileid)>1) THEN
                str30 = 't_' // tax_name(pfileid, itax)
             ELSE
                str30 = 'time_counter'
             END IF
             assoc = trim(str30) // ' ' // trim(assoc)

             iret = nf90_put_att(ncid, ncvarid, 'interval_operation', &
                  real(freq_opp(pfileid, iv)))
             iret = nf90_put_att(ncid, ncvarid, 'interval_write', real(freq_wrt( &
                  pfileid, iv)))
          END IF
          iret = nf90_put_att(ncid, ncvarid, 'associate', trim(assoc))
       END IF
    END DO

    !  Add MPP attributes

    CALL ioipslmpp_addatt(ncid)

    ! 3.0 Put the netcdf file into wrte mode

    iret = nf90_enddef(ncid)

    ! 4.0 Give some informations to the user

    WRITE (str70, '("All variables have been initialized on file:", I3)') &
         pfileid
    CALL histerr(1, 'histend', str70, '', ' ')

  END SUBROUTINE histend

  !**********************************

  SUBROUTINE histsync(file)

    ! This subroutine will synchronise all
    ! (or one if defined) opened files.

    ! file: optional argument for fileid

    USE histcom_var, ONLY: nb_files, ncdf_ids
    USE netcdf, ONLY: nf90_sync

    INTEGER, INTENT (IN), OPTIONAL:: file

    INTEGER:: ifile, ncid, iret

    LOGICAL:: file_exists
    !---------------------------------------------------------------------

    ! 1.The loop on files to synchronise

    DO ifile = 1, nb_files

       IF (present(file)) THEN
          file_exists = (ifile==file)
       ELSE
          file_exists = .TRUE.
       END IF

       IF (file_exists) THEN
          ncid = ncdf_ids(ifile)
          iret = nf90_sync(ncid)
       END IF

    END DO

  END SUBROUTINE histsync

  !**********************************

  SUBROUTINE histclo(fid)

    ! This subroutine will close the file corresponding
    ! to the argument, if the argument is present. Else it will close
    ! all opened files.

    USE errioipsl, ONLY: histerr
    USE histcom_var, ONLY: nb_files, ncdf_ids
    USE netcdf, ONLY: nf90_close, nf90_noerr

    INTEGER, INTENT (IN), OPTIONAL:: fid ! file id

    ! Variables local to the procedure:
    INTEGER ifile, ncid, iret, iv, ncvarid
    INTEGER start_loop, end_loop
    CHARACTER(len=70) str70

    !---------------------------------------------------------------------

    IF (present(fid)) THEN
       start_loop = fid
       end_loop = fid
    ELSE
       start_loop = 1
       end_loop = nb_files
    END IF

    DO ifile = start_loop, end_loop
       ncid = ncdf_ids(ifile)
       iret = nf90_close(ncid)
       IF (iret/=nf90_noerr) THEN
          WRITE(str70, '("This file has already been closed:", I3)') ifile
          CALL histerr(2, 'histclo', str70, '', '')
       END IF
    END DO

  END SUBROUTINE histclo

END MODULE histcom
