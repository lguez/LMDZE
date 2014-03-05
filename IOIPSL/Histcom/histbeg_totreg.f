MODULE histbeg_totreg_m

  ! From histcom.f90, version 2.1 2004/04/21 09:27:10

  ! Some confusing vocabulary in this code!

  ! A regular grid is a grid with i, j indices and thus it is
  ! stored in a 2D matrix. This is opposed to an irregular grid which
  ! is in a vector and where we do not know which neighbours we
  ! have. As a consequence we need the bounds for each grid-cell.

  ! A rectilinear grid is a special case of a regular grid in which
  ! all longitudes for i constant are equal and all latitudes for j
  ! constant are equal. In other words we do not need the full 2D
  ! matrix to describe the grid, just two vectors.

  IMPLICIT NONE

CONTAINS

  SUBROUTINE histbeg_totreg(filename, lon_1d, lat_1d, orix, szx, oriy, szy, &
       pitau0, pdate0, pdeltat, horiid, fileid)

    ! We assume the grid is rectilinear. The user provides "lon_1d"
    ! and "lat_1d" as vectors. This subroutine initializes a NetCDF
    ! file and returns the ID. It sets up the geographical space on
    ! which the data will be stored and offers the possibility of
    ! setting a zoom. It also gets the global parameters into the
    ! input-output subsystem.

    USE ioipslmpp, ONLY: ioipslmpp_file
    USE errioipsl, ONLY: histerr
    USE histcom_var, ONLY: assc_file, date0, deltat, full_size, itau0, &
         lock_modname, model_name, nb_files, nb_files_max, nb_hax, nb_tax, &
         nb_var, nb_zax, ncdf_ids, regular, slab_ori, slab_sz, xid, yid, zoom
    use histhori_regular_m, only: histhori_regular
    USE netcdf, ONLY: nf90_clobber, nf90_global
    use netcdf95, only: nf95_create, nf95_def_dim, nf95_put_att

    CHARACTER(len=*), INTENT(IN):: filename
    ! name of the netcdf file to be created

    REAL, INTENT(IN):: lon_1d(:) ! coordinates of points in longitude
    REAL, INTENT(IN):: lat_1d(:) ! coordinates of points in latitude

    ! The next 4 arguments allow to define a horizontal zoom for this
    ! file. It is assumed that all variables to come have the same
    ! index space. This can not be assumed for the z axis and thus we
    ! define the zoom in histdef.
    INTEGER, INTENT(IN):: orix ! origin of the slab of data within the X axis
    INTEGER, INTENT(IN):: szx ! size of the slab of data in X
    INTEGER, INTENT(IN):: oriy ! origin of the slab of data within the Y axis
    INTEGER, INTENT(IN):: szy ! size of the slab of data in Y

    INTEGER, INTENT(IN):: pitau0 ! time step at which the history tape starts
    REAL, INTENT(IN):: pdate0 ! the Julian date at which the itau was equal to 0
    REAL, INTENT(IN):: pdeltat ! time step of the counter itau, in seconds

    INTEGER, INTENT(OUT):: fileid ! ID of the netcdf file
    INTEGER, INTENT(OUT):: horiid ! ID of the horizontal grid

    ! Variables local to the procedure:
    REAL, DIMENSION(size(lon_1d), size(lat_1d)):: lon, lat
    INTEGER im ! size of arrays in longitude direction
    integer jm ! size of arrays in latitude direction
    INTEGER ncid
    INTEGER lengf, lenga
    CHARACTER(len=120) file

    !---------------------------------------------------------------------

    im = size(lon_1d)
    jm = size(lat_1d)

    lon = spread(lon_1d, 2, jm)
    lat = spread(lat_1d, 1, im)

    nb_files = nb_files + 1
    fileid = nb_files

    ! 1. Transfer into module variables for future use

    itau0(fileid) = pitau0
    date0(fileid) = pdate0
    deltat(fileid) = pdeltat

    ! 2. Initialize all variables for this file

    IF (nb_files > nb_files_max) CALL histerr(3, 'histbeg', &
         'Table of files too small. You should increase nb_files_max', &
         'in M_HISTCOM.f90 in order to accomodate all these files', ' ')

    nb_var(fileid) = 0
    nb_tax(fileid) = 0
    nb_hax(fileid) = 0
    nb_zax(fileid) = 0

    slab_ori(fileid, :) = (/orix, oriy/)
    slab_sz(fileid, :) = (/szx, szy/)

    ! 3. Open NetCDF file and define dimensions

    lengf = len_trim(filename)
    IF (filename(lengf-2:lengf)/='.nc') THEN
       file = filename(:lengf) // '.nc'
    ELSE
       file = filename(:lengf)
    END IF

    ! Add PE number in file name on MPP

    CALL ioipslmpp_file(file)

    ! Keep track of the name of the files opened

    lengf = len_trim(file)
    lenga = len_trim(assc_file)
    IF (nb_files==1) THEN
       assc_file = file(:lengf)
    ELSE IF ((lenga+lengf)<500) THEN
       assc_file = assc_file(:lenga) // ' ' // file(:lengf)
    ELSE IF (lenga + 7 < 500 .AND. index(assc_file(:lenga), 'et.al.') < 1) THEN
       assc_file = assc_file(:lenga) // ' et.al.'
    ELSE
       CALL histerr(2, 'histbeg', &
            'The file names do not fit into the associate_file attribute.', &
            'Use shorter names if you wish to keep the information.', ' ')
    END IF

    call nf95_create(file, nf90_clobber, ncid)
    call nf95_def_dim(ncid, 'lon', szx, xid(nb_files))
    call nf95_def_dim(ncid, 'lat', szy, yid(nb_files))

    ! 4. Declare the geographical coordinates and other attributes

    ! 4.3 Global attributes

    call nf95_put_att(ncid, nf90_global, 'Conventions', 'GDT 1.3')
    call nf95_put_att(ncid, nf90_global, 'file_name', trim(file))
    call nf95_put_att(ncid, nf90_global, 'production', trim(model_name))
    lock_modname = .TRUE.

    ! 5. Save some important information on this file in the module variables
    ncdf_ids(fileid) = ncid
    full_size(fileid, :) = (/im, jm/)

    ! 6. Store the geographical coordinates

    IF ((im /= szx) .OR. (jm /= szy)) zoom(fileid) = .TRUE.
    regular(fileid) = .TRUE.

    CALL histhori_regular(fileid, im, lon, jm, lat, ' ', 'Default grid', horiid)

  END SUBROUTINE histbeg_totreg

end MODULE histbeg_totreg_m
