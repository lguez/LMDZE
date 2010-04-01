MODULE flincom

  ! From flincom.f90, version 2.2 2006/03/07 09:21:51

  USE netcdf

  USE calendar,  ONLY : ju2ymds, ymds2ju, ioconf_calendar
  USE errioipsl, ONLY : histerr
  USE stringop,  ONLY : strlowercase

  IMPLICIT NONE

  PRIVATE
  PUBLIC flinput, flincre, flinget, flinget_zoom2d, flinclo
  public flinopen_nozoom
  public flininfo, flininspect, flinquery_var

  INTERFACE flinput
     !---------------------------------------------------------------------
     !- The "flinput" routines will put a variable
     !- on the netCDF file created by flincre.
     !- If the sizes of the axis do not match the one of the IDs
     !- then a new axis is created.
     !- That is we loose the possibility of writting hyperslabs of data.

     !- Again here if iim = jjm = llm = ttm = 0
     !- then a global attribute is added to the file.

     !- INPUT

     !- fid      : Identification of the file in which we will write
     !- varname  : Name of variable to be written
     !- iim      : size in x of variable
     !- nlonid   : ID of x axis which could fit for this axis
     !- jjm      : size in y of variable
     !- nlatid   : ID of y axis which could fit for this axis
     !- llm      : size in z of variable
     !- zdimid   : ID of z axis which could fit for this axis
     !- ttm      : size in t of variable
     !- tdimid   : ID of t axis which could fit for this axis

     !- OUTPUT

     !- NONE
     !---------------------------------------------------------------------
     MODULE PROCEDURE flinput_r4d, flinput_r3d, flinput_r2d, &
          flinput_r1d, flinput_scal
  END INTERFACE

  INTERFACE flinget
     MODULE PROCEDURE flinget_r4d, flinget_r3d, flinget_r2d, flinget_r1d, &
          flinget_scal
  END INTERFACE
  INTERFACE flinget_zoom2d
     MODULE PROCEDURE flinget_r4d_zoom2d, flinget_r3d_zoom2d, &
          flinget_r2d_zoom2d
  END INTERFACE

  ! This is the data we keep on each file we open

  INTEGER, PARAMETER :: nbfile_max = 200
  INTEGER, SAVE :: nbfiles = 0
  INTEGER, SAVE :: ncids(nbfile_max), ncnbd(nbfile_max), &
       ncfunli(nbfile_max), ncnba(nbfile_max)
  INTEGER, SAVE :: ncnbva(nbfile_max), ncdims(nbfile_max,4)
  LOGICAL, SAVE :: ncfileopen(nbfile_max)=.FALSE.

  INTEGER, SAVE :: cind_vid, cind_fid, cind_len
  INTEGER,DIMENSION(:),ALLOCATABLE,SAVE :: cindex

  INTEGER,DIMENSION(4) :: w_sta, w_len, w_dim

CONTAINS

  SUBROUTINE flincre &
       (filename, iim1, jjm1, lon1, lat1, llm1, lev1, ttm1, itaus, &
       time0, dt, fid_out, nlonid1, nlatid1, zdimid1, tdimid1)
    !---------------------------------------------------------------------
    !- This is a "low level" subroutine for opening netCDF files wich
    !- contain the major coordinate system of the model.
    !- Other coordinates needed for other variables
    !- will be added as they are needed.

    !- INPUT

    !- filename    : Name of the file to be created
    !- iim1, jjm1  : Horizontal size of the grid
    !-               which will be stored in the file
    !- lon1, lat1  : Horizontal grids
    !- llm1        : Size of the vertical grid
    !- lev1        : Vertical grid
    !- ttm1        : Size of time axis
    !- itaus       : time steps on the time axis
    !- time0       : Time in julian days at which itau = 0
    !- dt          : time step in seconds between itaus
    !-               (one step of itau)

    !- OUTPUT

    !- fid         : File identification
    !- nlonid1     : Identification of longitudinal axis
    !- nlatid1     : Identification of latitudinal axis
    !- zdimid1     : ID of vertical axis
    !- tdimid1     : ID of time axis
    !---------------------------------------------------------------------
    IMPLICIT NONE

    ! ARGUMENTS

    CHARACTER(LEN=*) :: filename
    INTEGER :: iim1, jjm1, llm1, ttm1
    REAL :: lon1(iim1,jjm1)
    REAL :: lat1(iim1,jjm1)
    REAL :: lev1(llm1)
    INTEGER :: itaus(ttm1)
    REAL :: time0
    REAL :: dt
    INTEGER :: fid_out, zdimid1, nlonid1, nlatid1, tdimid1

    ! LOCAL

    INTEGER :: iret, lll, fid
    INTEGER :: lonid, latid, levid, timeid
    INTEGER :: year, month, day
    REAL :: sec
    CHARACTER(LEN=250):: name

    LOGICAL :: check = .FALSE.
    !---------------------------------------------------------------------
    lll = LEN_TRIM(filename)
    IF (filename(lll-2:lll) /= '.nc') THEN
       name=filename(1:lll)//'.nc'
    ELSE
       name=filename(1:lll)
    ENDIF

    iret = NF90_CREATE (name, NF90_CLOBBER, fid)

    iret = NF90_DEF_DIM (fid, 'x',     iim1, nlonid1)
    iret = NF90_DEF_DIM (fid, 'y',     jjm1, nlatid1)
    iret = NF90_DEF_DIM (fid, 'lev',   llm1, zdimid1)
    iret = NF90_DEF_DIM (fid, 'tstep', ttm1, tdimid1)

    ! Vertical axis

    IF (check) WRITE(*,*) 'flincre Vertical axis'

    iret = NF90_DEF_VAR (fid, 'lev', NF90_FLOAT, zdimid1, levid)
    iret = NF90_PUT_ATT (fid, levid, 'units',     '-')
    iret = NF90_PUT_ATT (fid, levid, 'title',     'levels')
    iret = NF90_PUT_ATT (fid, levid, 'long_name', 'Sigma Levels')

    ! Time axis

    IF (check) WRITE(*,*) 'flincre time axis'

    iret = NF90_DEF_VAR (fid, 'tstep', NF90_FLOAT, tdimid1, timeid)
    iret = NF90_PUT_ATT (fid, timeid, 'units',     '-')
    iret = NF90_PUT_ATT (fid, timeid, 'title',     'time')
    iret = NF90_PUT_ATT (fid, timeid, 'long_name', 'time steps')

    ! The longitude

    IF (check) WRITE(*,*) 'flincre Longitude axis'

    iret = NF90_DEF_VAR (fid, "nav_lon", NF90_FLOAT, &
         (/ nlonid1, nlatid1 /), lonid)
    iret = NF90_PUT_ATT (fid, lonid, 'units', "degrees_east")
    iret = NF90_PUT_ATT (fid, lonid, 'title', "Longitude")
    iret = NF90_PUT_ATT (fid, lonid, 'nav_model', &
         "Lambert projection of PROMES")
    iret = NF90_PUT_ATT (fid, lonid, 'valid_min', &
         REAL(MINVAL(lon1)))
    iret = NF90_PUT_ATT (fid, lonid, 'valid_max', &
         REAL(MAXVAL(lon1)))

    ! The Latitude

    IF (check) WRITE(*,*) 'flincre Latitude axis'

    iret = NF90_DEF_VAR (fid, "nav_lat", NF90_FLOAT, &
         (/ nlonid1, nlatid1 /), latid)
    iret = NF90_PUT_ATT (fid, latid, 'units', "degrees_north")
    iret = NF90_PUT_ATT (fid, latid, 'title', "Latitude")
    iret = NF90_PUT_ATT (fid, latid, 'nav_model', &
         "Lambert projection of PROMES")
    iret = NF90_PUT_ATT (fid, latid, 'valid_min', &
         REAL(MINVAL(lat1)))
    iret = NF90_PUT_ATT (fid, latid, 'valid_max', &
         REAL(MAXVAL(lat1)))

    ! The time coordinates

    iret = NF90_PUT_ATT (fid, NF90_GLOBAL, 'delta_tstep_sec', &
         REAL(dt))

    CALL ju2ymds (time0, year, month, day, sec)

    iret = NF90_PUT_ATT (fid, NF90_GLOBAL, 'year0',  REAL(year))
    iret = NF90_PUT_ATT (fid, NF90_GLOBAL, 'month0', REAL(month))
    iret = NF90_PUT_ATT (fid, NF90_GLOBAL, 'day0',   REAL(day))
    iret = NF90_PUT_ATT (fid, NF90_GLOBAL, 'sec0',   REAL(sec))

    iret = NF90_ENDDEF (fid)

    IF (check) WRITE(*,*) 'flincre Variable'

    iret = NF90_PUT_VAR (fid, levid, lev1(1:llm1))

    IF (check) WRITE(*,*) 'flincre Time Variable'

    iret = NF90_PUT_VAR (fid, timeid, REAL(itaus(1:ttm1)))

    IF (check) WRITE(*,*) 'flincre Longitude'

    iret = NF90_PUT_VAR (fid, lonid, lon1(1:iim1,1:jjm1))

    IF (check) WRITE(*,*) 'flincre Latitude'

    iret = NF90_PUT_VAR (fid, latid, lat1(1:iim1,1:jjm1))

    ! Keep all this information

    nbfiles = nbfiles+1

    IF (nbfiles > nbfile_max) THEN
       CALL histerr (3,'flincre', &
            'Too many files. Please increase nbfil_max', &
            'in program flincom.F90.',' ')
    ENDIF

    ncids(nbfiles) = fid
    ncnbd(nbfiles) = 4

    ncdims(nbfiles,1:4) = (/ iim1, jjm1, llm1, ttm1 /)

    ncfunli(nbfiles) = -1
    ncnba(nbfiles)   =  4
    ncnbva(nbfiles)  =  0
    ncfileopen(nbfiles) = .TRUE.

    fid_out = nbfiles
    !---------------------
  END SUBROUTINE flincre

  !===

  SUBROUTINE flinopen_nozoom(filename, iim, jjm, llm, lon, lat, lev, &
       ttm, itaus, date0, dt, fid_out)

    !- The routine will open an input file
    !- INPUT
    !- filename  : Name of the netCDF file to be opened

    !- There is no test of the content of the file against the input
    ! from the model

    !- iim       : size in the x direction in the file (longitude)
    !- jjm       : size in the y direction
    !- llm       : number of levels
    !-             (llm = 0 means no axis to be expected)

    !- WARNING :
    !- It is for the user to check
    !- that the dimensions of lon lat and lev are correct when passed to
    !- flinopen. This can be done after the call when iim and jjm have
    !- been retrieved from the netCDF file. In F90 this problem will
    !- be solved with an internal assign
    !- IF iim, jjm, llm or ttm are parameters in the calling program 
    !- it will create a segmentation fault

    !-  ttm       : size of time axis

    !-  OUTPUT

    !- lon       : array of (iim,jjm),
    !-             that contains the longitude of each point
    !- lat       : same for latitude
    !- lev       : An array of llm for the latitude
    !-  itaus     : Time steps within this file
    !-  date0     : Julian date at which itau = 0
    !-  dt        : length of the time steps of the data
    
    !---------------------------------------------------------------------

    IMPLICIT NONE

    ! ARGUMENTS

    CHARACTER(LEN=*), intent(in):: filename
    INTEGER, intent(in) :: iim, jjm, llm, ttm
    real, intent(out):: lon(iim,jjm), lat(iim,jjm), lev(llm)
    INTEGER, intent(out):: itaus(ttm)
    REAL, intent(out):: date0, dt

    INTEGER, intent(in):: fid_out
    ! (file ID which is later used to read the data)

    ! LOCAL

    INTEGER :: iret, vid, fid, nbdim, i
    INTEGER :: gdtt_id, old_id, iv, gdtmaf_id
    CHARACTER(LEN=250) :: name
    CHARACTER(LEN=80) :: units, calendar
    INTEGER :: year, month, day
    REAL :: r_year, r_month, r_day
    INTEGER :: year0, month0, day0, hours0, minutes0, seci
    REAL :: sec, sec0
    CHARACTER :: strc

    REAL,DIMENSION(:),ALLOCATABLE :: vec_tmp

    !---------------------------------------------------------------------

    IF ( (fid_out < 1).OR.(fid_out > nbfile_max) ) THEN
       ! Either the fid_out has not been initialized (0 or very large)
       ! then we have to open anyway. Else we only need to open the file
       ! if it has not been opened before.
       print *, "Call flinfo before flinopen"
       stop 1
    end IF
    IF (.NOT.ncfileopen(fid_out)) THEN
       print *, "Call flinfo before flinopen"
       stop 1
    end IF

    !-- The user has already opened the file
    !-- and we trust that he knows the dimensions

    fid = ncids(fid_out)

    ! 2.0 get the sizes and names of the different coordinates
    !     and do a first set of verification.

    ! 3.0 Check if we are realy talking about the same coodinate system
    !     if not then we get the lon, lat and lev variables from the file

    !-- 4.0 extracting the coordinates
    !---
    CALL flinfindcood (fid_out, 'lon', vid, nbdim)
    IF (nbdim == 2) THEN
       iret = NF90_GET_VAR (fid, vid, lon, &
            start=(/ 1, 1 /), count=(/ iim, jjm /))
    ELSE
       ALLOCATE(vec_tmp(iim))
       iret = NF90_GET_VAR (fid, vid, vec_tmp, &
            start=(/ 1 /), count=(/ iim /))
       DO i=1,jjm
          lon(:,i) = vec_tmp(:)
       ENDDO
       DEALLOCATE(vec_tmp)
    ENDIF
    !---
    CALL flinfindcood (fid_out, 'lat', vid, nbdim)
    IF (nbdim == 2) THEN
       iret = NF90_GET_VAR (fid, vid, lat, &
            start=(/ 1, 1 /), count=(/ iim, jjm /))
    ELSE
       ALLOCATE(vec_tmp(jjm))
       iret = NF90_GET_VAR (fid, vid, vec_tmp, &
            start=(/ 1 /), count=(/ jjm /))
       DO i=1,iim
          lat(i,:) = vec_tmp(:)
       ENDDO
       DEALLOCATE(vec_tmp)
    ENDIF
    !---
    IF (llm > 0) THEN
       CALL flinfindcood (fid_out, 'lev', vid, nbdim)
       IF (nbdim == 1) THEN
          iret = NF90_GET_VAR (fid, vid, lev, &
               start=(/ 1 /), count=(/ llm /))
       ELSE
          CALL histerr (3,'flinopen', &
               'Can not handle vertical coordinates that have more',&
               'than 1 dimension',' ')
       ENDIF
    ENDIF

    ! 5.0 Get all the details for the time if possible needed

    IF (ttm > 0) THEN
       !---
       !-- 5.1 Find the time axis. Prefered method is the 'timestep since'
       !---
       gdtmaf_id = -1
       gdtt_id = -1
       old_id = -1
       DO iv=1,ncnbva(fid_out)
          name=''
          iret = NF90_INQUIRE_VARIABLE (fid, iv, name=name)
          units=''
          iret = NF90_GET_ATT (fid, iv, 'units', units)
          IF (INDEX(units,'seconds since') > 0) gdtmaf_id = iv
          IF (INDEX(units,'timesteps since') > 0) gdtt_id = iv
          IF (INDEX(name, 'tstep') > 0) old_id = iv
       ENDDO
       !---
       IF (gdtt_id > 0) THEN
          vid = gdtt_id
       ELSE IF (gdtmaf_id > 0) THEN
          vid = gdtmaf_id
       ELSE IF (old_id > 0) THEN
          vid = old_id
       ELSE
          CALL histerr (3, 'flinopen', 'No time axis found',' ',' ')
       ENDIF
       !---
       ALLOCATE(vec_tmp(ttm))
       iret = NF90_GET_VAR (fid, vid, vec_tmp, &
            start=(/ 1 /), count=(/ ttm /))
       itaus(1:ttm) = NINT(vec_tmp(1:ttm))
       DEALLOCATE(vec_tmp)
       !---
       !-- Getting all the details for the time axis
       !---
       !-- Find the calendar
       calendar='XXXX'
       iret = NF90_GET_ATT (fid, gdtmaf_id, 'calendar', calendar)
       IF ( INDEX(calendar,'XXXX') < 1 ) THEN
          CALL ioconf_calendar(calendar)
       ENDIF
       !--
       units = ''
       iret = NF90_GET_ATT (fid, vid, 'units', units)
       IF (gdtt_id > 0) THEN
          units = units(INDEX(units,'since')+6:LEN_TRIM(units))
          READ (units,'(I4.4,5(a,I2.2))') &
               year0, strc, month0, strc, day0, &
               strc, hours0, strc, minutes0, strc, seci
          sec0 = hours0*3600. + minutes0*60. + seci
          CALL ymds2ju (year0, month0, day0, sec0, date0)
          iret = NF90_GET_ATT (fid, gdtt_id, 'tstep_sec', dt)
       ELSE IF (gdtmaf_id > 0) THEN
          units = units(INDEX(units,'since')+6:LEN_TRIM(units))
          READ (units,'(I4.4,5(a,I2.2))') &
               year0, strc, month0, strc, day0, &
               strc, hours0, strc, minutes0, strc, seci
          sec0 = hours0*3600. + minutes0*60. + seci
          CALL ymds2ju (year0, month0, day0, sec0, date0)
          !-----
       ELSE IF (old_id > 0) THEN
          iret = NF90_GET_ATT (fid, NF90_GLOBAL, 'delta_tstep_sec', dt)
          iret = NF90_GET_ATT (fid, NF90_GLOBAL, 'day0', r_day)
          iret = NF90_GET_ATT (fid, NF90_GLOBAL, 'sec0', sec)
          iret = NF90_GET_ATT (fid, NF90_GLOBAL, 'year0', r_year)
          iret = NF90_GET_ATT (fid, NF90_GLOBAL, 'month0', r_month)
          !-----
          day = INT(r_day)
          month = INT(r_month)
          year = INT(r_year)
          !-----
          CALL ymds2ju (year, month, day, sec, date0)
       ENDIF
    ENDIF

  END SUBROUTINE flinopen_nozoom

  !===

  SUBROUTINE flininfo (filename, iim, jjm, llm, ttm, fid_out)
    !---------------------------------------------------------------------
    !- This subroutine allows to get some information.
    !- It is usualy done within flinopen but the user may want to call
    !- it before in order to allocate the space needed to extract the
    !- data from the file.
    !---------------------------------------------------------------------
    IMPLICIT NONE

    ! ARGUMENTS

    CHARACTER(LEN=*), intent(in):: filename
    INTEGER, intent(out):: iim, jjm, llm, ttm, fid_out

    ! LOCAL

    INTEGER :: iret, fid, ndims, nvars, nb_atts, id_unlim
    INTEGER :: iv, lll
    CHARACTER(LEN=80) :: name
    CHARACTER(LEN=30) :: axname

    LOGICAL :: check = .FALSE.
    !---------------------------------------------------------------------
    lll = LEN_TRIM(filename)
    IF (filename(lll-2:lll) /= '.nc') THEN
       name = filename(1:lll)//'.nc'
    ELSE
       name = filename(1:lll)
    ENDIF

    iret = NF90_OPEN (name, NF90_NOWRITE, fid)
    IF (iret /= NF90_NOERR) THEN
       CALL histerr(3, 'flininfo','Could not open file :',TRIM(name),' ')
    ENDIF

    iret = NF90_INQUIRE (fid, nDimensions=ndims, nVariables=nvars, &
         nAttributes=nb_atts, unlimitedDimId=id_unlim)

    iim = 0;
    jjm = 0;
    llm = 0;
    ttm = 0;

    DO iv=1,ndims
       !---
       iret = NF90_INQUIRE_DIMENSION (fid, iv, name=axname, len=lll)
       CALL strlowercase (axname)
       axname = ADJUSTL(axname)
       !---
       IF (check) WRITE(*,*) &
            'flininfo - getting axname',iv,axname,lll
       !---
       IF      (    (INDEX(axname,'x') == 1) &
            .OR.(INDEX(axname,'lon') == 1) ) THEN
          iim = lll;
       ELSE IF (    (INDEX(axname,'y') == 1) &
            .OR.(INDEX(axname,'lat') == 1) ) THEN
          jjm = lll;
       ELSE IF (    (INDEX(axname,'lev') == 1) &
            .OR.(INDEX(axname,'plev') == 1) &
            .OR.(INDEX(axname,'z') == 1) &
            .OR.(INDEX(axname,'depth') == 1) ) THEN
          llm = lll;
       ELSE IF (    (INDEX(axname,'tstep') == 1) &
            .OR.(INDEX(axname,'time_counter') == 1) ) THEN
          !---- For the time we certainly need to allow for other names
          ttm = lll;
       ELSE IF (ndims == 1) THEN
          !---- Nothing was found and ndims=1 then we have a vector of data
          iim = lll;
       ENDIF
       !---
    ENDDO

    ! Keep all this information

    nbfiles = nbfiles+1

    IF (nbfiles > nbfile_max) THEN
       CALL histerr (3,'flininfo', &
            'Too many files. Please increase nbfil_max', &
            'in program flincom.F90.',' ')
    ENDIF

    ncids(nbfiles) = fid
    ncnbd(nbfiles) = ndims

    ncdims(nbfiles,1:4) = (/ iim, jjm, llm, ttm /)

    ncfunli(nbfiles) = id_unlim
    ncnba(nbfiles)   = nb_atts
    ncnbva(nbfiles)  = nvars
    ncfileopen(nbfiles) = .TRUE.

    fid_out = nbfiles
    !----------------------
  END SUBROUTINE flininfo

  !===

  SUBROUTINE flinput_r1d &
       (fid_in,varname,iim,nlonid,jjm,nlatid,llm,zdimid,ttm,tdimid,var)
    !---------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER :: fid_in
    CHARACTER(LEN=*) :: varname
    INTEGER :: iim, nlonid, jjm, nlatid, llm, zdimid, ttm, tdimid
    REAL :: var(:)

    INTEGER :: fid, ncvarid, ndim, iret
    LOGICAL :: check = .FALSE.
    !---------------------------------------------------------------------
    IF (check) WRITE(*,*) &
         "flinput_r1d : SIZE(var) = ",SIZE(var)

    CALL flinput_mat &
         (fid_in,varname,iim,nlonid,jjm,nlatid,llm,zdimid,ttm,tdimid, &
         fid,ncvarid,ndim)

    iret = NF90_PUT_VAR (fid, ncvarid, var, &
         start=w_sta(1:ndim), count=w_len(1:ndim))
    !-------------------------
  END SUBROUTINE flinput_r1d

  !===

  SUBROUTINE flinput_r2d &
       (fid_in,varname,iim,nlonid,jjm,nlatid,llm,zdimid,ttm,tdimid,var)
    !---------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER :: fid_in
    CHARACTER(LEN=*) :: varname
    INTEGER :: iim, nlonid, jjm, nlatid, llm, zdimid, ttm, tdimid
    REAL :: var(:,:)

    INTEGER :: fid, ncvarid, ndim, iret
    LOGICAL :: check = .FALSE.
    !---------------------------------------------------------------------
    IF (check) WRITE(*,*) &
         "flinput_r2d : SIZE(var) = ",SIZE(var)

    CALL flinput_mat &
         (fid_in,varname,iim,nlonid,jjm,nlatid,llm,zdimid,ttm,tdimid, &
         fid,ncvarid,ndim)

    iret = NF90_PUT_VAR (fid, ncvarid, var, &
         start=w_sta(1:ndim), count=w_len(1:ndim))
    !-------------------------
  END SUBROUTINE flinput_r2d

  !===

  SUBROUTINE flinput_r3d &
       (fid_in,varname,iim,nlonid,jjm,nlatid,llm,zdimid,ttm,tdimid,var)
    !---------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER :: fid_in
    CHARACTER(LEN=*) :: varname
    INTEGER :: iim, nlonid, jjm, nlatid, llm, zdimid, ttm, tdimid
    REAL :: var(:,:,:)

    INTEGER :: fid, ncvarid, ndim, iret
    LOGICAL :: check = .FALSE.
    !---------------------------------------------------------------------
    IF (check) WRITE(*,*) &
         "flinput_r3d : SIZE(var) = ",SIZE(var)

    CALL flinput_mat &
         (fid_in,varname,iim,nlonid,jjm,nlatid,llm,zdimid,ttm,tdimid, &
         fid,ncvarid,ndim)

    iret = NF90_PUT_VAR (fid, ncvarid, var, &
         start=w_sta(1:ndim), count=w_len(1:ndim))
    !-------------------------
  END SUBROUTINE flinput_r3d

  !===

  SUBROUTINE flinput_r4d &
       (fid_in,varname,iim,nlonid,jjm,nlatid,llm,zdimid,ttm,tdimid,var)
    !---------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER :: fid_in
    CHARACTER(LEN=*) :: varname
    INTEGER :: iim, nlonid, jjm, nlatid, llm, zdimid, ttm, tdimid
    REAL :: var(:,:,:,:)

    INTEGER :: fid, ncvarid, ndim, iret
    LOGICAL :: check = .FALSE.
    !---------------------------------------------------------------------
    IF (check) WRITE(*,*) &
         "flinput_r4d : SIZE(var) = ",SIZE(var)

    CALL flinput_mat &
         (fid_in,varname,iim,nlonid,jjm,nlatid,llm,zdimid,ttm,tdimid, &
         fid,ncvarid,ndim)

    iret = NF90_PUT_VAR (fid, ncvarid, var, &
         start=w_sta(1:ndim), count=w_len(1:ndim))
    !-------------------------
  END SUBROUTINE flinput_r4d

  !===

  SUBROUTINE flinput_mat &
       (fid_in,varname,iim,nlonid,jjm,nlatid, &
       llm,zdimid,ttm,tdimid,fid,ncvarid,ndim)
    !---------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER :: fid_in
    CHARACTER(LEN=*) :: varname
    INTEGER :: iim, nlonid, jjm, nlatid, llm, zdimid, ttm, tdimid
    INTEGER :: fid, ncvarid, ndim

    ! LOCAL

    INTEGER :: iret
    !---------------------------------------------------------------------
    fid = ncids(fid_in)

    w_sta(1:4) = (/      1,      1,  1,  1 /)
    w_len(1:2) = (/    iim,    jjm /)
    w_dim(1:2) = (/ nlonid, nlatid /)

    IF ( (llm > 0).AND.(ttm > 0) ) THEN
       ndim = 4
       w_len(3:4) = (/    llm,    ttm /)
       w_dim(3:4) = (/ zdimid, tdimid /)
    ELSE IF (llm > 0) THEN
       ndim = 3
       w_dim(3) = zdimid
       w_len(3) = llm
    ELSE IF (ttm > 0) THEN
       ndim = 3
       w_dim(3) = tdimid
       w_len(3) = ttm
    ELSE
       ndim = 2
    ENDIF

    iret = NF90_REDEF   (fid)
    iret = NF90_DEF_VAR (fid,varname,NF90_FLOAT,w_dim(1:ndim),ncvarid)
    iret = NF90_PUT_ATT (fid,ncvarid,'short_name',TRIM(varname))
    iret = NF90_ENDDEF  (fid)
    !--------------------------
  END  SUBROUTINE flinput_mat

  !===

  SUBROUTINE flinput_scal &
       (fid_in, varname, iim, nlonid, jjm, nlatid, &
       llm, zdimid, ttm, tdimid, var)
    !---------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER :: fid_in
    CHARACTER(LEN=*) :: varname
    INTEGER :: iim, nlonid, jjm, nlatid, llm, zdimid, ttm, tdimid
    REAL :: var

    ! LOCAL

    INTEGER :: fid, iret
    !---------------------------------------------------------------------
    fid = ncids(fid_in)

    iret = NF90_REDEF   (fid)
    iret = NF90_PUT_ATT (fid, NF90_GLOBAL, varname, REAL(var))
    iret = NF90_ENDDEF  (fid)
    !---------------------------
  END  SUBROUTINE flinput_scal

  !===

  SUBROUTINE flinget_r1d &
       (fid_in,varname,iim,jjm,llm,ttm,itau_dep,itau_fin,var)
    !---------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER :: fid_in
    CHARACTER(LEN=*) :: varname
    INTEGER :: iim, jjm, llm, ttm, itau_dep, itau_fin
    REAL :: var(:)

    INTEGER :: jl, ji
    REAL,DIMENSION(:),ALLOCATABLE,SAVE :: buff_tmp
    LOGICAL :: check = .FALSE.
    !---------------------------------------------------------------------
    IF (.NOT.ALLOCATED(buff_tmp)) THEN
       IF (check) WRITE(*,*) &
            "flinget_r1d : allocate buff_tmp for buff_sz = ",SIZE(var)
       ALLOCATE (buff_tmp(SIZE(var)))
    ELSE IF (SIZE(var) > SIZE(buff_tmp)) THEN
       IF (check) WRITE(*,*) &
            "flinget_r1d : re-allocate buff_tmp for buff_sz = ",SIZE(var)
       DEALLOCATE (buff_tmp)
       ALLOCATE (buff_tmp(SIZE(var)))
    ENDIF

    CALL flinget_mat (fid_in,varname,iim,jjm,llm,ttm, &
         itau_dep,itau_fin,1,iim,1,jjm,buff_tmp)

    jl=0
    DO ji=1,SIZE(var,1)
       jl=jl+1
       var(ji) = buff_tmp(jl)
    ENDDO
    !-------------------------
  END SUBROUTINE flinget_r1d

  !===

  SUBROUTINE flinget_r2d &
       (fid_in,varname,iim,jjm,llm,ttm,itau_dep,itau_fin,var)
    !---------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER :: fid_in
    CHARACTER(LEN=*) :: varname
    INTEGER :: iim, jjm, llm, ttm, itau_dep, itau_fin
    REAL :: var(:,:)

    INTEGER :: jl, jj, ji
    REAL,DIMENSION(:),ALLOCATABLE,SAVE :: buff_tmp
    LOGICAL :: check = .FALSE.
    !---------------------------------------------------------------------
    IF (.NOT.ALLOCATED(buff_tmp)) THEN
       IF (check) WRITE(*,*) &
            "flinget_r2d : allocate buff_tmp for buff_sz = ",SIZE(var)
       ALLOCATE (buff_tmp(SIZE(var)))
    ELSE IF (SIZE(var) > SIZE(buff_tmp)) THEN
       IF (check) WRITE(*,*) &
            "flinget_r2d : re-allocate buff_tmp for buff_sz = ",SIZE(var)
       DEALLOCATE (buff_tmp)
       ALLOCATE (buff_tmp(SIZE(var)))
    ENDIF

    CALL flinget_mat (fid_in,varname,iim,jjm,llm,ttm, &
         itau_dep,itau_fin,1,iim,1,jjm,buff_tmp)

    jl=0
    DO jj=1,SIZE(var,2)
       DO ji=1,SIZE(var,1)
          jl=jl+1
          var(ji,jj) = buff_tmp(jl)
       ENDDO
    ENDDO
    !-------------------------
  END SUBROUTINE flinget_r2d

  !===

  SUBROUTINE flinget_r2d_zoom2d &
       (fid_in,varname,iim,jjm,llm,ttm, &
       itau_dep,itau_fin,iideb,iilen,jjdeb,jjlen,var)
    !---------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER :: fid_in
    CHARACTER(LEN=*) :: varname
    INTEGER :: iim,jjm,llm,ttm,itau_dep,itau_fin,iideb,jjdeb,iilen,jjlen
    REAL :: var(:,:)

    INTEGER :: jl, jj, ji
    REAL,DIMENSION(:),ALLOCATABLE,SAVE :: buff_tmp
    LOGICAL :: check = .FALSE.
    !---------------------------------------------------------------------
    IF (.NOT.ALLOCATED(buff_tmp)) THEN
       IF (check) WRITE(*,*) &
            "flinget_r2d_zoom : allocate buff_tmp for buff_sz = ",SIZE(var)
       ALLOCATE (buff_tmp(SIZE(var)))
    ELSE IF (SIZE(var) > SIZE(buff_tmp)) THEN
       IF (check) WRITE(*,*) &
            "flinget_r2d_zoom : re-allocate buff_tmp for buff_sz = ",SIZE(var)
       DEALLOCATE (buff_tmp)
       ALLOCATE (buff_tmp(SIZE(var)))
    ENDIF

    CALL flinget_mat (fid_in,varname,iim,jjm,llm,ttm, &
         itau_dep,itau_fin,iideb,iilen,jjdeb,jjlen,buff_tmp)

    jl=0
    DO jj=1,SIZE(var,2)
       DO ji=1,SIZE(var,1)
          jl=jl+1
          var(ji,jj) = buff_tmp(jl)
       ENDDO
    ENDDO
    !--------------------------------
  END SUBROUTINE flinget_r2d_zoom2d

  !===

  SUBROUTINE flinget_r3d(fid_in,varname,iim,jjm,llm,ttm,itau_dep,itau_fin,var)
    !---------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER, intent(in):: fid_in
    CHARACTER(LEN=*), intent(in):: varname
    INTEGER, intent(in):: iim, jjm, llm, ttm, itau_dep, itau_fin
    REAL, intent(out):: var(:,:,:)

    INTEGER :: jl, jk, jj, ji
    REAL,DIMENSION(:),ALLOCATABLE,SAVE :: buff_tmp
    LOGICAL :: check = .FALSE.
    !---------------------------------------------------------------------
    IF (.NOT.ALLOCATED(buff_tmp)) THEN
       IF (check) WRITE(*,*) &
            "flinget_r3d : allocate buff_tmp for buff_sz = ",SIZE(var)
       ALLOCATE (buff_tmp(SIZE(var)))
    ELSE IF (SIZE(var) > SIZE(buff_tmp)) THEN
       IF (check) WRITE(*,*) &
            "flinget_r3d : re-allocate buff_tmp for buff_sz = ",SIZE(var)
       DEALLOCATE (buff_tmp)
       ALLOCATE (buff_tmp(SIZE(var)))
    ENDIF

    CALL flinget_mat (fid_in,varname,iim,jjm,llm,ttm, &
         itau_dep,itau_fin,1,iim,1,jjm,buff_tmp)

    jl=0
    DO jk=1,SIZE(var,3)
       DO jj=1,SIZE(var,2)
          DO ji=1,SIZE(var,1)
             jl=jl+1
             var(ji,jj,jk) = buff_tmp(jl)
          ENDDO
       ENDDO
    ENDDO
    !-------------------------
  END SUBROUTINE flinget_r3d

  !===

  SUBROUTINE flinget_r3d_zoom2d &
       (fid_in,varname,iim,jjm,llm,ttm, &
       itau_dep,itau_fin,iideb,iilen,jjdeb,jjlen,var)
    !---------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER :: fid_in
    CHARACTER(LEN=*) :: varname
    INTEGER :: iim,jjm,llm,ttm,itau_dep,itau_fin,iideb,jjdeb,iilen,jjlen
    REAL :: var(:,:,:)

    INTEGER :: jl, jk, jj, ji
    REAL,DIMENSION(:),ALLOCATABLE,SAVE :: buff_tmp
    LOGICAL :: check = .FALSE.
    !---------------------------------------------------------------------
    IF (.NOT.ALLOCATED(buff_tmp)) THEN
       IF (check) WRITE(*,*) &
            "flinget_r3d_zoom : allocate buff_tmp for buff_sz = ",SIZE(var)
       ALLOCATE (buff_tmp(SIZE(var)))
    ELSE IF (SIZE(var) > SIZE(buff_tmp)) THEN
       IF (check) WRITE(*,*) &
            "flinget_r3d_zoom : re-allocate buff_tmp for buff_sz = ",SIZE(var)
       DEALLOCATE (buff_tmp)
       ALLOCATE (buff_tmp(SIZE(var)))
    ENDIF

    CALL flinget_mat (fid_in,varname,iim,jjm,llm,ttm, &
         itau_dep,itau_fin,iideb,iilen,jjdeb,jjlen,buff_tmp)

    jl=0
    DO jk=1,SIZE(var,3)
       DO jj=1,SIZE(var,2)
          DO ji=1,SIZE(var,1)
             jl=jl+1
             var(ji,jj,jk) = buff_tmp(jl)
          ENDDO
       ENDDO
    ENDDO
    !--------------------------------
  END SUBROUTINE flinget_r3d_zoom2d

  !===

  SUBROUTINE flinget_r4d &
       (fid_in,varname,iim,jjm,llm,ttm,itau_dep,itau_fin,var)
    !---------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER :: fid_in
    CHARACTER(LEN=*) :: varname
    INTEGER :: iim, jjm, llm, ttm, itau_dep, itau_fin
    REAL :: var(:,:,:,:)

    INTEGER :: jl, jk, jj, ji, jm
    REAL,DIMENSION(:),ALLOCATABLE,SAVE :: buff_tmp
    LOGICAL :: check = .FALSE.
    !---------------------------------------------------------------------
    IF (.NOT.ALLOCATED(buff_tmp)) THEN
       IF (check) WRITE(*,*) &
            "flinget_r4d : allocate buff_tmp for buff_sz = ",SIZE(var)
       ALLOCATE (buff_tmp(SIZE(var)))
    ELSE IF (SIZE(var) > SIZE(buff_tmp)) THEN
       IF (check) WRITE(*,*) &
            "flinget_r4d : re-allocate buff_tmp for buff_sz = ",SIZE(var)
       DEALLOCATE (buff_tmp)
       ALLOCATE (buff_tmp(SIZE(var)))
    ENDIF

    CALL flinget_mat (fid_in,varname,iim,jjm,llm,ttm, &
         itau_dep,itau_fin,1,iim,1,jjm,buff_tmp)

    jl=0
    DO jm=1,SIZE(var,4)
       DO jk=1,SIZE(var,3)
          DO jj=1,SIZE(var,2)
             DO ji=1,SIZE(var,1)
                jl=jl+1
                var(ji,jj,jk,jm) = buff_tmp(jl)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !-------------------------
  END SUBROUTINE flinget_r4d

  !===

  SUBROUTINE flinget_r4d_zoom2d &
       (fid_in,varname,iim,jjm,llm,ttm, &
       itau_dep,itau_fin,iideb,iilen,jjdeb,jjlen,var)
    !---------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER :: fid_in
    CHARACTER(LEN=*) :: varname
    INTEGER :: iim,jjm,llm,ttm,itau_dep,itau_fin,iideb,jjdeb,iilen,jjlen
    REAL :: var(:,:,:,:)

    INTEGER :: jl, jk, jj, ji, jm
    REAL,DIMENSION(:),ALLOCATABLE,SAVE :: buff_tmp
    LOGICAL :: check = .FALSE.
    !---------------------------------------------------------------------
    IF (.NOT.ALLOCATED(buff_tmp)) THEN
       IF (check) WRITE(*,*) &
            "flinget_r4d_zoom : allocate buff_tmp for buff_sz = ",SIZE(var)
       ALLOCATE (buff_tmp(SIZE(var)))
    ELSE IF (SIZE(var) > SIZE(buff_tmp)) THEN
       IF (check) WRITE(*,*) &
            "flinget_r4d_zoom : re-allocate buff_tmp for buff_sz = ",SIZE(var)
       DEALLOCATE (buff_tmp)
       ALLOCATE (buff_tmp(SIZE(var)))
    ENDIF

    CALL flinget_mat (fid_in,varname,iim,jjm,llm,ttm, &
         itau_dep,itau_fin,iideb,iilen,jjdeb,jjlen,buff_tmp)

    jl=0
    DO jm=1,SIZE(var,4)
       DO jk=1,SIZE(var,3)
          DO jj=1,SIZE(var,2)
             DO ji=1,SIZE(var,1)
                jl=jl+1
                var(ji,jj,jk,jm) = buff_tmp(jl)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !--------------------------------
  END SUBROUTINE flinget_r4d_zoom2d

  !===

  SUBROUTINE flinget_mat &
       (fid_in, varname, iim, jjm, llm, ttm, itau_dep, &
       itau_fin, iideb, iilen, jjdeb, jjlen, var)
    !---------------------------------------------------------------------
    !- This subroutine will read the variable named varname from
    !- the file previously opened by flinopen and identified by fid

    !- It is checked that the dimensions of the variable to be read
    !- correspond to what the user requested when he specified
    !- iim, jjm and llm. The only exception which is allowed is
    !- for compressed data where the horizontal grid is not expected
    !- to be iim x jjm.

    !- If variable is of size zero a global attribute is read.
    !- This global attribute will be of type real

    !- INPUT

    !- fid      : File ID returned by flinopen
    !- varname  : Name of the variable to be read from the file
    !- iim      : | These three variables give the size of the variables
    !- jjm      : | to be read. It will be verified that the variables
    !- llm      : | fits in there.
    !- ttm      : |
    !- itau_dep : Time step at which we will start to read
    !- itau_fin : Time step until which we are going to read
    !-            For the moment this is done on indexes
    !-            but it should be in the physical space.
    !-            If there is no time-axis in the file then use a
    !-            itau_fin < itau_dep, this will tell flinget not to
    !-            expect a time-axis in the file.
    !- iideb    : index i for zoom
    !- iilen    : length of zoom
    !- jjdeb    : index j for zoom
    !- jjlen    : length of zoom

    !- OUTPUT

    !- var      : array that will contain the data
    !---------------------------------------------------------------------
    IMPLICIT NONE

    ! ARGUMENTS

    INTEGER, intent(in):: fid_in
    CHARACTER(LEN=*), intent(in):: varname
    INTEGER, intent(in):: iim, jjm, llm, ttm, itau_dep, itau_fin
    INTEGER :: iideb
    integer, intent(in):: iilen
    integer jjdeb
    integer, intent(in):: jjlen
    REAL :: var(:)

    ! LOCAL

    INTEGER :: iret, fid
    INTEGER :: vid, cvid, clen
    CHARACTER(LEN=70) :: str1
    CHARACTER(LEN=250) :: att_n, tmp_n
    CHARACTER(LEN=5) :: axs_l
    INTEGER :: tmp_i
    REAL,SAVE :: mis_v=0.
    REAL :: tmp_r
    INTEGER :: ndims, x_typ, nb_atts
    INTEGER,DIMENSION(NF90_MAX_VAR_DIMS) :: dimids
    INTEGER :: i, iv, nvars, i2d, cnd
    REAL,DIMENSION(:),ALLOCATABLE,SAVE :: var_tmp
    LOGICAL :: uncompress = .FALSE.
    LOGICAL :: check = .FALSE.
    !---------------------------------------------------------------------
    fid = ncids(fid_in)

    IF (check) THEN
       WRITE(*,*) &
            'flinget_mat : fid_in, fid, varname :', fid_in, fid, TRIM(varname)
       WRITE(*,*) &
            'flinget_mat : iim, jjm, llm, ttm, itau_dep, itau_fin :', &
            iim, jjm, llm, ttm, itau_dep, itau_fin
       WRITE(*,*) &
            'flinget_mat : iideb, iilen, jjdeb, jjlen :', &
            iideb, iilen, jjdeb, jjlen
    ENDIF

    uncompress = .FALSE.

    ! 1.0 We get first all the details on this variable from the file

    nvars = ncnbva(fid_in)

    vid = -1
    iret = NF90_INQ_VARID (fid, varname, vid)

    IF (vid < 0 .OR. iret /= NF90_NOERR) THEN
       CALL histerr (3,'flinget', &
            'Variable '//TRIM(varname)//' not found in file',' ',' ')
    ENDIF

    iret = NF90_INQUIRE_VARIABLE (fid, vid, &
         ndims=ndims, dimids=dimids, nAtts=nb_atts)
    IF (check) THEN
       WRITE(*,*) &
            'flinget_mat : fid, vid :', fid, vid
       WRITE(*,*) &
            'flinget_mat : ndims, dimids(1:ndims), nb_atts :', &
            ndims, dimids(1:ndims), nb_atts
    ENDIF

    w_dim(:) = 0
    DO i=1,ndims
       iret  = NF90_INQUIRE_DIMENSION (fid, dimids(i), len=w_dim(i))
    ENDDO
    IF (check) WRITE(*,*) &
         'flinget_mat : w_dim :', w_dim(1:ndims)

    mis_v = 0.0; axs_l = ' ';

    IF (nb_atts > 0) THEN
       IF (check) THEN
          WRITE(*,*) 'flinget_mat : attributes for variable :'
       ENDIF
    ENDIF
    DO i=1,nb_atts
       iret = NF90_INQ_ATTNAME (fid, vid, i, att_n)
       iret = NF90_INQUIRE_ATTRIBUTE (fid, vid, att_n, xtype=x_typ)
       CALL strlowercase (att_n)
       IF      (    (x_typ == NF90_INT).OR.(x_typ == NF90_SHORT) &
            .OR.(x_typ == NF90_BYTE) ) THEN
          iret = NF90_GET_ATT (fid, vid, att_n, tmp_i)
          IF (check) THEN
             WRITE(*,*) '   ',TRIM(att_n),' : ',tmp_i
          ENDIF
       ELSE IF ( (x_typ == NF90_FLOAT).OR.(x_typ == NF90_DOUBLE) ) THEN
          iret = NF90_GET_ATT (fid, vid, att_n, tmp_r)
          IF (check) THEN
             WRITE(*,*) '   ',TRIM(att_n),' : ',tmp_r
          ENDIF
          IF (index(att_n,'missing_value') > 0) THEN
             mis_v = tmp_r
          ENDIF
       ELSE
          tmp_n = ''
          iret = NF90_GET_ATT (fid, vid, att_n, tmp_n)
          IF (check) THEN
             WRITE(*,*) '   ',TRIM(att_n),' : ',TRIM(tmp_n)
          ENDIF
          IF (index(att_n,'axis') > 0) THEN
             axs_l = tmp_n
          ENDIF
       ENDIF
    ENDDO
    !?
!!!!!!!!!! We will need a verification on the type of the variable
    !?

    ! 2.0 The dimensions are analysed to determine what is to be read

    ! 2.1 the longitudes

    IF ( w_dim(1) /= iim .OR. w_dim(2) /= jjm) THEN
       !---
       !-- There is a possibility that we have to deal with a compressed axis !
       !---
       iret = NF90_INQUIRE_DIMENSION (fid, dimids(1), &
            name=tmp_n, len=clen)
       iret = NF90_INQ_VARID (fid, tmp_n, cvid)
       !---
       IF (check) WRITE(*,*) &
            'Dimname, iret , NF90_NOERR : ',TRIM(tmp_n),iret,NF90_NOERR
       !---
       !-- If we have an axis which has the same name
       !-- as the dimension we can see if it is compressed
       !---
       !-- TODO TODO for zoom2d
       !---
       IF (iret == NF90_NOERR) THEN
          iret = NF90_GET_ATT (fid, cvid, 'compress', str1)
          !-----
          IF (iret == NF90_NOERR) THEN
             iret = NF90_INQUIRE_VARIABLE (fid,cvid,xtype=x_typ,ndims=cnd)
             !-------
             IF ( cnd /= 1 .AND. x_typ /= NF90_INT) THEN
                CALL histerr (3,'flinget', &
                     'Variable '//TRIM(tmp_n)//' can not be a compressed axis', &
                     'Either it has too many dimensions'// &
                     ' or it is not of type integer', ' ')
             ELSE
                !---------
                !-------- Let us see if we already have that index table
                !---------
                IF (    (cind_len /= clen).OR.(cind_vid /= cvid) &
                     .OR.(cind_fid /= fid) ) THEN
                   IF (ALLOCATED(cindex))   DEALLOCATE(cindex)
                   ALLOCATE(cindex(clen))
                   cind_len = clen
                   cind_vid = cvid
                   cind_fid = fid
                   iret = NF90_GET_VAR (fid, cvid, cindex)
                ENDIF
                !---------
                !-------- In any case we need to set the slab of data to be read
                !---------
                uncompress = .TRUE.
                w_sta(1) = 1
                w_len(1) = clen
                i2d = 1
             ENDIF
          ELSE
             str1 = 'The horizontal dimensions of '//varname
             CALL histerr (3,'flinget',str1, &
                  'is not compressed and does not'// &
                  ' correspond to the requested size',' ')
          ENDIF
       ELSE
          IF (w_dim(1) /= iim) THEN
             str1 = 'The longitude dimension of '//varname
             CALL histerr (3,'flinget',str1, &
                  'in the file is not equal to the dimension', &
                  'that should be read')
          ENDIF
          IF (w_dim(2) /= jjm) THEN
             str1 = 'The latitude dimension of '//varname
             CALL histerr (3,'flinget',str1, &
                  'in the file is not equal to the dimension', &
                  'that should be read')
          ENDIF
       ENDIF
    ELSE
       w_sta(1:2) = (/ iideb, jjdeb /)
       w_len(1:2) = (/ iilen, jjlen /)
       i2d = 2
    ENDIF

    ! 2.3 Now the difficult part, the 3rd dimension which can be
    ! time or levels.

    ! Priority is given to the time axis if only three axes are present.

    IF (ndims > i2d) THEN
       !---
       !-- 2.3.1 We have a vertical axis
       !---
       IF (llm == 1 .AND. ndims == i2d+2 .OR. llm == w_dim(i2d+1)) THEN
          !-----
          IF (w_dim(i2d+1) /= llm) THEN
             CALL histerr (3,'flinget', &
                  'The vertical dimension of '//varname, &
                  'in the file is not equal to the dimension', &
                  'that should be read')
          ELSE
             w_sta(i2d+1) = 1
             IF (llm > 0) THEN
                w_len(i2d+1) = llm
             ELSE
                w_len(i2d+1) = w_sta(i2d+1)
             ENDIF
          ENDIF
          !-----
          IF ((itau_fin-itau_dep) >= 0) THEN
             IF      (ndims /= i2d+2) THEN
                CALL histerr (3,'flinget', &
                     'You attempt to read a time slab', &
                     'but there is no time axis on this variable', varname)
             ELSE IF ((itau_fin - itau_dep) <= w_dim(i2d+2)) THEN
                w_sta(i2d+2) = itau_dep
                w_len(i2d+2) = itau_fin-itau_dep+1
             ELSE
                CALL histerr (3,'flinget', &
                     'The time step you try to read is not', &
                     'in the file (1)', varname)
             ENDIF
          ELSE IF (ndims == i2d+2 .AND. w_dim(i2d+2) > 1) THEN
             CALL histerr (3,'flinget', &
                  'There is a time axis in the file but no', &
                  'time step give in the call', varname)
          ELSE
             w_sta(i2d+2) = 1
             w_len(i2d+2) = 1
          ENDIF
       ELSE
          !-----
          !---- 2.3.2 We do not have any vertical axis
          !-----
          IF (ndims == i2d+2) THEN
             CALL histerr (3,'flinget', &
                  'The file contains 4 dimensions', &
                  'but only 3 are requestes for variable ', varname)
          ENDIF
          IF ((itau_fin-itau_dep) >= 0) THEN
             IF (ndims == i2d+1) THEN
                IF ((itau_fin-itau_dep) < w_dim(i2d+1) ) THEN
                   w_sta(i2d+1) = itau_dep
                   w_len(i2d+1) = itau_fin-itau_dep+1
                ELSE
                   CALL histerr (3,'flinget', &
                        'The time step you try to read is not', &
                        'in the file (2)', varname)
                ENDIF
             ELSE
                CALL histerr (3,'flinget', &
                     'From your input you sould have 3 dimensions', &
                     'in the file but there are 4', varname)
             ENDIF
          ELSE
             IF (ndims == i2d+1 .AND. w_dim(i2d+1) > 1) THEN
                CALL histerr (3,'flinget', &
                     'There is a time axis in the file but no', &
                     'time step given in the call', varname)
             ELSE
                w_sta(i2d+1) = 1
                w_len(i2d+1) = 1
             ENDIF
          ENDIF
       ENDIF
    ELSE
       !---
       !-- 2.3.3 We do not have any vertical axis
       !---
       w_sta(i2d+1:i2d+2) = (/ 0, 0 /)
       w_len(i2d+1:i2d+2) = (/ 0, 0 /)
    ENDIF

    ! 3.0 Reading the data

    IF (check) WRITE(*,*) &
         'flinget_mat 3.0 : ', uncompress, w_sta, w_len
    !---
    IF (uncompress) THEN
       !---
       IF (ALLOCATED(var_tmp)) THEN
          IF (SIZE(var_tmp) < clen) THEN
             DEALLOCATE(var_tmp)
             ALLOCATE(var_tmp(clen))
          ENDIF
       ELSE
          ALLOCATE(var_tmp(clen))
       ENDIF
       !---
       iret = NF90_GET_VAR (fid, vid, var_tmp, &
            start=w_sta(:), count=w_len(:))
       !---
       var(:) = mis_v
       var(cindex(:)) = var_tmp(:)
       !---
    ELSE
       iret = NF90_GET_VAR (fid, vid, var, &
            start=w_sta(:), count=w_len(:))
    ENDIF

    IF (check) WRITE(*,*) 'flinget_mat 3.1 : ',NF90_STRERROR (iret)
    !--------------------------
  END  SUBROUTINE flinget_mat

  !===

  SUBROUTINE flinget_scal &
       (fid_in, varname, iim, jjm, llm, ttm, itau_dep, itau_fin, var)
    !---------------------------------------------------------------------
    !- This subroutine will read the variable named varname from
    !- the file previously opened by flinopen and identified by fid

    !- If variable is of size zero a global attribute is read. This
    !- global attribute will be of type real

    !- INPUT

    !- fid      : File ID returned by flinopen
    !- varname  : Name of the variable to be read from the file
    !- iim      : | These three variables give the size of the variables
    !- jjm      : | to be read. It will be verified that the variables
    !- llm      : | fits in there.
    !- ttm      : |
    !- itau_dep : Time step at which we will start to read
    !- itau_fin : Time step until which we are going to read
    !-           For the moment this is done on indeces but it should be
    !-           in the physical space
    !-           If there is no time-axis in the file then use a
    !-           itau_fin < itau_dep, this will tell flinget not to
    !-           expect a time-axis in the file.

    !- OUTPUT

    !- var      : scalar that will contain the data
    !---------------------------------------------------------------------
    IMPLICIT NONE

    ! ARGUMENTS

    INTEGER :: fid_in
    CHARACTER(LEN=*) :: varname
    INTEGER :: iim, jjm, llm, ttm, itau_dep, itau_fin
    REAL :: var

    ! LOCAL

    INTEGER :: iret, fid

    LOGICAL :: check = .FALSE.
    !---------------------------------------------------------------------
    fid = ncids(fid_in)

    ! 1.0 Reading a global attribute

    iret = NF90_GET_ATT (fid, NF90_GLOBAL, varname, var)
    !---------------------------
  END  SUBROUTINE flinget_scal

  !===

  SUBROUTINE flinfindcood (fid_in, axtype, vid, ndim)
    !---------------------------------------------------------------------
    !- This subroutine explores the file in order to find
    !- the coordinate according to a number of rules
    !---------------------------------------------------------------------
    IMPLICIT NONE

    ! ARGUMENTS

    INTEGER, intent(in):: fid_in
    integer vid, ndim
    CHARACTER(LEN=3) :: axtype

    ! LOCAL

    INTEGER :: iv, iret, dimnb
    CHARACTER(LEN=40) :: dimname, dimuni1, dimuni2, dimuni3
    CHARACTER(LEN=30) :: str1
    LOGICAL :: found_rule = .FALSE.
    !---------------------------------------------------------------------
    vid = -1

    ! Make sure all strings are invalid

    dimname = '?-?'
    dimuni1 = '?-?'
    dimuni2 = '?-?'
    dimuni3 = '?-?'

    ! First rule : we look for the correct units
    ! lon : east
    ! lat : north
    ! We make an exact check as it would be too easy to mistake
    ! some units by just comparing the substrings.

    SELECTCASE(axtype)
    CASE ('lon')
       dimuni1 = 'degree_e'
       dimuni2 = 'degrees_e'
       found_rule = .TRUE.
    CASE('lat')
       dimuni1 = 'degree_n'
       dimuni2 = 'degrees_n'
       found_rule = .TRUE.
    CASE('lev')
       dimuni1 = 'm'
       dimuni2 = 'km'
       dimuni3 = 'hpa'
       found_rule = .TRUE.
    CASE DEFAULT
       found_rule = .FALSE.
    END SELECT

    IF (found_rule) THEN
       iv = 0
       DO WHILE ( (vid < 0).AND.(iv < ncnbva(fid_in)) )
          iv = iv+1
          str1 = ''
          iret = NF90_GET_ATT (ncids(fid_in), iv, 'units', str1)
          IF (iret == NF90_NOERR) THEN
             CALL strlowercase (str1)
             IF (    (INDEX(str1, TRIM(dimuni1)) == 1) &
                  .OR.(INDEX(str1, TRIM(dimuni2)) == 1) &
                  .OR.(INDEX(str1, TRIM(dimuni3)) == 1) ) THEN
                vid = iv
                iret = NF90_INQUIRE_VARIABLE (ncids(fid_in), iv, ndims=ndim)
             ENDIF
          ENDIF
       ENDDO
    ENDIF

    ! Second rule : we find specific names :
    ! lon : nav_lon
    ! lat : nav_lat
    ! Here we can check if we find the substring as the
    ! names are more specific.

    SELECTCASE(axtype)
    CASE ('lon')
       dimname = 'nav_lon lon longitude'
       found_rule = .TRUE.
    CASE('lat')
       dimname = 'nav_lat lat latitude'
       found_rule = .TRUE.
    CASE('lev')
       dimname = 'plev level depth deptht'
       found_rule = .TRUE.
    CASE DEFAULT
       found_rule = .FALSE.
    END SELECT

    IF (found_rule) THEN
       iv = 0
       DO WHILE ( (vid < 0).AND.(iv < ncnbva(fid_in)) )
          iv = iv+1
          str1=''
          iret = NF90_INQUIRE_VARIABLE (ncids(fid_in), iv, &
               name=str1, ndims=ndim)
          IF (INDEX(dimname,TRIM(str1)) >= 1) THEN
             vid = iv
          ENDIF
       ENDDO
    ENDIF

    ! Third rule : we find a variable with the same name as the dimension
    ! lon = 1
    ! lat = 2
    ! lev = 3

    IF (vid < 0) THEN
       SELECTCASE(axtype)
       CASE ('lon')
          dimnb = 1
          found_rule = .TRUE.
       CASE('lat')
          dimnb = 2
          found_rule = .TRUE.
       CASE('lev')
          dimnb = 3
          found_rule = .TRUE.
       CASE DEFAULT
          found_rule = .FALSE.
       END SELECT
       !---
       IF (found_rule) THEN
          iret = NF90_INQUIRE_DIMENSION (ncids(fid_in), dimnb, name=dimname)
          iv = 0
          DO WHILE ( (vid < 0).AND.(iv < ncnbva(fid_in)) )
             iv = iv+1
             str1=''
             iret = NF90_INQUIRE_VARIABLE (ncids(fid_in), iv, &
                  name=str1, ndims=ndim)
             IF (INDEX(dimname,TRIM(str1)) == 1) THEN
                vid = iv
             ENDIF
          ENDDO
       ENDIF
    ENDIF

    ! Stop the program if no coordinate was found

    IF (vid < 0) THEN
       CALL histerr (3,'flinfindcood', &
            'No coordinate axis was found in the file', &
            'The data in this file can not be used', axtype)
    ENDIF
    !--------------------------
  END SUBROUTINE flinfindcood

  !===

  SUBROUTINE flinclo (fid_in)
    !---------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER :: fid_in

    INTEGER :: iret
    !---------------------------------------------------------------------
    iret = NF90_CLOSE (ncids(fid_in))
    ncfileopen(fid_in) = .FALSE.
    !---------------------
  END SUBROUTINE flinclo

  !===

  SUBROUTINE flinquery_var(fid_in, varname, exists)
    !---------------------------------------------------------------------
    !- Queries the existance of a variable in the file.
    !---------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER :: fid_in
    CHARACTER(LEN=*) varname
    LOGICAL :: exists

    INTEGER :: iret, fid, vid
    !---------------------------------------------------------------------
    fid = ncids(fid_in)
    vid = -1
    iret = NF90_INQ_VARID (fid, varname, vid)

    exists = ( (vid >= 0).AND.(iret == NF90_NOERR) )
    !---------------------------
  END SUBROUTINE flinquery_var

  !===

  SUBROUTINE flininspect (fid_in)
    !---------------------------------------------------------------------
    IMPLICIT NONE

    ! fid : File id to inspect

    INTEGER :: fid_in

    !- LOCAL

    INTEGER :: iim, jjm, llm, ttm
    INTEGER :: iret, fid, ndims, nvars, nb_atts, id_unlim
    INTEGER :: iv, in, lll
    INTEGER :: xid, yid, zid, tid
    INTEGER,DIMENSION(NF90_MAX_VAR_DIMS) :: idimid
    CHARACTER(LEN=80) :: name
    CHARACTER(LEN=30) :: axname
    !---------------------------------------------------------------------
    fid = ncids(fid_in)

    iret = NF90_INQUIRE (fid, nDimensions=ndims, nVariables=nvars, &
         nAttributes=nb_atts, unlimitedDimId=id_unlim)

    WRITE (*,*) 'IOIPSL ID                   : ',fid_in
    WRITE (*,*) 'NetCDF ID                   : ',fid
    WRITE (*,*) 'Number of dimensions        : ',ndims
    WRITE (*,*) 'Number of variables         : ',nvars
    WRITE (*,*) 'Number of global attributes : ',nb_atts
    WRITE (*,*) 'ID unlimited                : ',id_unlim

    xid = -1; iim = 0;
    yid = -1; jjm = 0;
    zid = -1; llm = 0;
    tid = -1; ttm = 0;

    DO iv=1,ndims
       !---
       iret = NF90_INQUIRE_DIMENSION (fid, iv, name=axname, len=lll)
       CALL strlowercase (axname)
       axname = ADJUSTL(axname)
       !---
       WRITE (*,*) 'Dimension number : ',iv
       WRITE (*,*) 'Dimension name   : ',TRIM(axname)
       !---
       IF      (    (INDEX(axname,'x') == 1) &
            .OR.(INDEX(axname,'lon') == 1)) THEN
          xid = iv; iim = lll;
          WRITE (*,*) 'Dimension X size   : ',iim
       ELSE IF (    (INDEX(axname,'y') == 1) &
            .OR.(INDEX(axname,'lat') == 1)) THEN
          yid = iv; jjm = lll;
          WRITE (*,*) 'Dimension Y size   : ',jjm
       ELSE IF (    (INDEX(axname,'lev') == 1) &
            .OR.(INDEX(axname,'plev') == 1) &
            .OR.(INDEX(axname,'z') == 1) &
            .OR.(INDEX(axname,'depth') == 1)) THEN
          zid = iv; llm = lll;
          WRITE (*,*) 'Dimension Z size   : ',llm
       ELSE IF (    (INDEX(axname,'tstep') == 1) &
            .OR.(INDEX(axname,'time_counter') == 1)) THEN
          !---- For the time we certainly need to allow for other names
          tid = iv; ttm = lll;
       ELSE IF (ndims == 1) THEN
          !---- Nothing was found and ndims=1 then we have a vector of data
          xid = 1; iim = lll;
       ENDIF
       !---
    ENDDO

    ! Keep all this information

    nbfiles = nbfiles+1

    IF (nbfiles > nbfile_max) THEN
       CALL histerr(3,'flininspect', &
            'Too many files. Please increase nbfil_max', &
            'in program flincom.F90.',' ')
    ENDIF

    ncids(nbfiles) = fid
    ncnbd(nbfiles) = ndims

    ncdims(nbfiles,1:4) = (/ iim, jjm, llm, ttm /)

    ncfunli(nbfiles) = id_unlim
    ncnba(nbfiles)   = nb_atts
    ncnbva(nbfiles)  = nvars
    ncfileopen(nbfiles) = .TRUE.

    DO in=1,nvars
       iret = NF90_INQUIRE_VARIABLE (fid, in, &
            name=name, ndims=ndims, dimids=idimid, nAtts=nb_atts)
       WRITE (*,*) 'Variable number  ------------ > ', in
       WRITE (*,*) 'Variable name        : ', TRIM(name)
       WRITE (*,*) 'Number of dimensions : ', ndims
       WRITE (*,*) 'Dimensions ID''s     : ', idimid(1:ndims)
       WRITE (*,*) 'Number of attributes : ', nb_atts
    ENDDO
    !-------------------------
  END SUBROUTINE flininspect

  !===

END MODULE flincom
