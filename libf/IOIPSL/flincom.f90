MODULE flincom

  ! From flincom.f90, version 2.2 2006/03/07 09:21:51

  IMPLICIT NONE

  PRIVATE
  PUBLIC flinclo, flinopen_nozoom, flininfo, ncids

  ! This is the data we keep on each file we open:
  INTEGER, PARAMETER:: nbfile_max = 200
  INTEGER, SAVE:: ncids(nbfile_max)
  INTEGER, SAVE:: ncnbva(nbfile_max)
  LOGICAL, SAVE:: ncfileopen(nbfile_max)=.FALSE.

CONTAINS

  SUBROUTINE flinopen_nozoom(iim, jjm, llm, lon, lat, lev, &
       ttm, itaus, date0, dt, fid_out)

    ! The routine will open an input file
    ! INPUT
    ! There is no test of the content of the file against the input
    ! from the model

    ! iim: size in the x direction in the file (longitude)
    ! jjm: size in the y direction
    ! llm: number of levels
    ! (llm = 0 means no axis to be expected)

    ! WARNING:
    ! It is for the user to check
    ! that the dimensions of lon lat and lev are correct when passed to
    ! flinopen. This can be done after the call when iim and jjm have
    ! been retrieved from the netCDF file. In F90 this problem will
    ! be solved with an internal assign
    ! IF iim, jjm, llm or ttm are parameters in the calling program 
    ! it will create a segmentation fault

    ! ttm: size of time axis

    ! OUTPUT

    ! lon: array of (iim, jjm),
    ! that contains the longitude of each point
    ! lat: same for latitude
    ! lev: An array of llm for the latitude
    ! itaus: Time steps within this file
    ! date0: Julian date at which itau = 0
    ! dt: length of the time steps of the data

    !---------------------------------------------------------------------

    USE calendar, ONLY: ymds2ju, ioconf_calendar
    USE errioipsl, ONLY: histerr
    USE netcdf, ONLY: nf90_get_att, nf90_get_var, nf90_global, &
         nf90_inquire_variable

    ! ARGUMENTS

    INTEGER, intent(in):: iim, jjm, llm, ttm
    real, intent(out):: lon(iim, jjm), lat(iim, jjm), lev(llm)
    INTEGER, intent(out):: itaus(ttm)
    REAL, intent(out):: date0, dt

    INTEGER, intent(in):: fid_out
    ! (file ID which is later used to read the data)

    ! LOCAL

    INTEGER:: iret, vid, fid, nbdim, i
    INTEGER:: gdtt_id, old_id, iv, gdtmaf_id
    CHARACTER(LEN=250):: name
    CHARACTER(LEN=80):: units, my_calendar
    INTEGER:: year, month, day
    REAL:: r_year, r_month, r_day
    INTEGER:: year0, month0, day0, hours0, minutes0, seci
    REAL:: sec, sec0
    CHARACTER:: strc

    REAL, DIMENSION(:), ALLOCATABLE:: vec_tmp

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

    ! The user has already opened the file
    ! and we trust that he knows the dimensions

    fid = ncids(fid_out)

    ! 2.0 get the sizes and names of the different coordinates
    ! and do a first set of verification.

    ! 3.0 Check if we are realy talking about the same coodinate system
    ! if not then we get the lon, lat and lev variables from the file

    ! 4.0 extracting the coordinates

    CALL flinfindcood (fid_out, 'lon', vid, nbdim)
    IF (nbdim == 2) THEN
       iret = NF90_GET_VAR (fid, vid, lon, &
            start=(/ 1, 1 /), count=(/ iim, jjm /))
    ELSE
       ALLOCATE(vec_tmp(iim))
       iret = NF90_GET_VAR (fid, vid, vec_tmp, &
            start=(/ 1 /), count=(/ iim /))
       DO i=1, jjm
          lon(:, i) = vec_tmp(:)
       ENDDO
       DEALLOCATE(vec_tmp)
    ENDIF

    CALL flinfindcood (fid_out, 'lat', vid, nbdim)
    IF (nbdim == 2) THEN
       iret = NF90_GET_VAR (fid, vid, lat, &
            start=(/ 1, 1 /), count=(/ iim, jjm /))
    ELSE
       ALLOCATE(vec_tmp(jjm))
       iret = NF90_GET_VAR (fid, vid, vec_tmp, &
            start=(/ 1 /), count=(/ jjm /))
       DO i=1, iim
          lat(i, :) = vec_tmp(:)
       ENDDO
       DEALLOCATE(vec_tmp)
    ENDIF

    IF (llm > 0) THEN
       CALL flinfindcood (fid_out, 'lev', vid, nbdim)
       IF (nbdim == 1) THEN
          iret = NF90_GET_VAR (fid, vid, lev, &
               start=(/ 1 /), count=(/ llm /))
       ELSE
          CALL histerr (3, 'flinopen', &
               'Can not handle vertical coordinates that have more', &
               'than 1 dimension', ' ')
       ENDIF
    ENDIF

    ! 5.0 Get all the details for the time if possible needed

    IF (ttm > 0) THEN

       ! 5.1 Find the time axis. Prefered method is the 'timestep since'

       gdtmaf_id = -1
       gdtt_id = -1
       old_id = -1
       DO iv=1, ncnbva(fid_out)
          name=''
          iret = NF90_INQUIRE_VARIABLE (fid, iv, name=name)
          units=''
          iret = NF90_GET_ATT (fid, iv, 'units', units)
          IF (INDEX(units, 'seconds since') > 0) gdtmaf_id = iv
          IF (INDEX(units, 'timesteps since') > 0) gdtt_id = iv
          IF (INDEX(name, 'tstep') > 0) old_id = iv
       ENDDO

       IF (gdtt_id > 0) THEN
          vid = gdtt_id
       ELSE IF (gdtmaf_id > 0) THEN
          vid = gdtmaf_id
       ELSE IF (old_id > 0) THEN
          vid = old_id
       ELSE
          CALL histerr (3, 'flinopen', 'No time axis found', ' ', ' ')
       ENDIF

       ALLOCATE(vec_tmp(ttm))
       iret = NF90_GET_VAR (fid, vid, vec_tmp, &
            start=(/ 1 /), count=(/ ttm /))
       itaus(1:ttm) = NINT(vec_tmp(1:ttm))
       DEALLOCATE(vec_tmp)

       ! Getting all the details for the time axis

       ! Find the calendar
       my_calendar='XXXX'
       iret = NF90_GET_ATT (fid, gdtmaf_id, 'calendar', my_calendar)
       IF ( INDEX(my_calendar, 'XXXX') < 1 ) THEN
          CALL ioconf_calendar(my_calendar)
       ENDIF

       units = ''
       iret = NF90_GET_ATT (fid, vid, 'units', units)
       IF (gdtt_id > 0) THEN
          units = units(INDEX(units, 'since')+6:LEN_TRIM(units))
          READ (units, '(I4.4, 5(a, I2.2))') &
               year0, strc, month0, strc, day0, &
               strc, hours0, strc, minutes0, strc, seci
          sec0 = hours0*3600. + minutes0*60. + seci
          CALL ymds2ju (year0, month0, day0, sec0, date0)
          iret = NF90_GET_ATT (fid, gdtt_id, 'tstep_sec', dt)
       ELSE IF (gdtmaf_id > 0) THEN
          units = units(INDEX(units, 'since')+6:LEN_TRIM(units))
          READ (units, '(I4.4, 5(a, I2.2))') &
               year0, strc, month0, strc, day0, &
               strc, hours0, strc, minutes0, strc, seci
          sec0 = hours0*3600. + minutes0*60. + seci
          CALL ymds2ju (year0, month0, day0, sec0, date0)

       ELSE IF (old_id > 0) THEN
          iret = NF90_GET_ATT (fid, NF90_GLOBAL, 'delta_tstep_sec', dt)
          iret = NF90_GET_ATT (fid, NF90_GLOBAL, 'day0', r_day)
          iret = NF90_GET_ATT (fid, NF90_GLOBAL, 'sec0', sec)
          iret = NF90_GET_ATT (fid, NF90_GLOBAL, 'year0', r_year)
          iret = NF90_GET_ATT (fid, NF90_GLOBAL, 'month0', r_month)

          day = INT(r_day)
          month = INT(r_month)
          year = INT(r_year)

          CALL ymds2ju (year, month, day, sec, date0)
       ENDIF
    ENDIF

  END SUBROUTINE flinopen_nozoom

  !***************************************************************

  SUBROUTINE flininfo(filename, iim, jjm, llm, ttm, fid_out)

    ! This subroutine allows to get some information.
    ! It is usualy done within flinopen but the user may want to call
    ! it before in order to allocate the space needed to extract the
    ! data from the file.

    USE strlowercase_m, ONLY: strlowercase
    USE errioipsl, ONLY: histerr
    USE netcdf, ONLY: nf90_inquire, nf90_inquire_dimension, nf90_noerr, &
         nf90_nowrite
    use netcdf95, only: nf95_open

    CHARACTER(LEN=*), intent(in):: filename
    INTEGER, intent(out):: iim, jjm, llm, ttm, fid_out

    ! LOCAL

    INTEGER, SAVE:: nbfiles = 0
    INTEGER, SAVE:: ncdims(nbfile_max, 4)
    INTEGER:: iret, fid, ndims, nvars, nb_atts, id_unlim
    INTEGER:: iv, lll
    CHARACTER(LEN=80):: name
    CHARACTER(LEN=30):: axname

    !---------------------------------------------------------------------

    lll = LEN_TRIM(filename)
    IF (filename(lll-2:lll) /= '.nc') THEN
       name = filename(1:lll)//'.nc'
    ELSE
       name = filename(1:lll)
    ENDIF

    call NF95_OPEN(name, NF90_NOWRITE, fid)
    iret = NF90_INQUIRE(fid, nDimensions=ndims, nVariables=nvars, &
         nAttributes=nb_atts, unlimitedDimId=id_unlim)

    iim = 0
    jjm = 0
    llm = 0
    ttm = 0

    DO iv=1, ndims
       iret = NF90_INQUIRE_DIMENSION(fid, iv, name=axname, len=lll)
       CALL strlowercase(axname)
       axname = ADJUSTL(axname)

       IF ((INDEX(axname, 'x') == 1) .OR. (INDEX(axname, 'lon') == 1)) THEN
          iim = lll
       ELSE IF ((INDEX(axname, 'y') == 1) &
            .OR. (INDEX(axname, 'lat') == 1)) THEN
          jjm = lll
       ELSE IF ((INDEX(axname, 'lev') == 1) .OR. (INDEX(axname, 'plev') == 1) &
            .OR. (INDEX(axname, 'z') == 1) &
            .OR. (INDEX(axname, 'depth') == 1)) THEN
          llm = lll
       ELSE IF ((INDEX(axname, 'tstep') == 1) &
            .OR. (INDEX(axname, 'time_counter') == 1)) THEN
          ! For the time we certainly need to allow for other names
          ttm = lll
       ELSE IF (ndims == 1) THEN
          ! Nothing was found and ndims=1 then we have a vector of data
          iim = lll
       ENDIF
    ENDDO

    ! Keep all this information

    nbfiles = nbfiles+1

    IF (nbfiles > nbfile_max) THEN
       CALL histerr(3, 'flininfo', &
            'Too many files. Please increase nbfil_max', &
            'in program flincom.F90.', ' ')
    ENDIF

    ncids(nbfiles) = fid
    ncdims(nbfiles, :) = (/ iim, jjm, llm, ttm /)
    ncnbva(nbfiles) = nvars
    ncfileopen(nbfiles) = .TRUE.
    fid_out = nbfiles

  END SUBROUTINE flininfo

  !***************************************************************

  SUBROUTINE flinfindcood (fid_in, axtype, vid, ndim)

    ! This subroutine explores the file in order to find
    ! the coordinate according to a number of rules

    USE strlowercase_m, ONLY: strlowercase
    USE errioipsl, ONLY: histerr
    USE netcdf, ONLY: nf90_get_att, nf90_inquire_dimension, &
         nf90_inquire_variable, nf90_noerr

    ! ARGUMENTS

    INTEGER, intent(in):: fid_in
    integer vid, ndim
    CHARACTER(LEN=3):: axtype

    ! LOCAL

    INTEGER:: iv, iret, dimnb
    CHARACTER(LEN=40):: dimname, dimuni1, dimuni2, dimuni3
    CHARACTER(LEN=30):: str1
    LOGICAL:: found_rule = .FALSE.
    !---------------------------------------------------------------------
    vid = -1

    ! Make sure all strings are invalid

    dimname = '?-?'
    dimuni1 = '?-?'
    dimuni2 = '?-?'
    dimuni3 = '?-?'

    ! First rule: we look for the correct units
    ! lon: east
    ! lat: north
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
             IF ( (INDEX(str1, TRIM(dimuni1)) == 1) &
                  .OR.(INDEX(str1, TRIM(dimuni2)) == 1) &
                  .OR.(INDEX(str1, TRIM(dimuni3)) == 1) ) THEN
                vid = iv
                iret = NF90_INQUIRE_VARIABLE (ncids(fid_in), iv, ndims=ndim)
             ENDIF
          ENDIF
       ENDDO
    ENDIF

    ! Second rule: we find specific names:
    ! lon: nav_lon
    ! lat: nav_lat
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
          IF (INDEX(dimname, TRIM(str1)) >= 1) THEN
             vid = iv
          ENDIF
       ENDDO
    ENDIF

    ! Third rule: we find a variable with the same name as the dimension
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

       IF (found_rule) THEN
          iret = NF90_INQUIRE_DIMENSION (ncids(fid_in), dimnb, name=dimname)
          iv = 0
          DO WHILE ( (vid < 0).AND.(iv < ncnbva(fid_in)) )
             iv = iv+1
             str1=''
             iret = NF90_INQUIRE_VARIABLE (ncids(fid_in), iv, &
                  name=str1, ndims=ndim)
             IF (INDEX(dimname, TRIM(str1)) == 1) THEN
                vid = iv
             ENDIF
          ENDDO
       ENDIF
    ENDIF

    ! Stop the program if no coordinate was found

    IF (vid < 0) THEN
       CALL histerr (3, 'flinfindcood', &
            'No coordinate axis was found in the file', &
            'The data in this file can not be used', axtype)
    ENDIF

  END SUBROUTINE flinfindcood

  !***************************************************************

  SUBROUTINE flinclo (fid_in)

    USE netcdf, ONLY: nf90_close

    INTEGER:: fid_in

    INTEGER:: iret

    !---------------------------------------------------------------------

    iret = NF90_CLOSE (ncids(fid_in))
    ncfileopen(fid_in) = .FALSE.

  END SUBROUTINE flinclo

END MODULE flincom
