MODULE histbeg_totreg_m

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
    use histhori_regular_m, only: histhori_regular
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

end MODULE histbeg_totreg_m
