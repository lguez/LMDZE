module flinopen_nozoom_m

  implicit none

contains

  SUBROUTINE flinopen_nozoom(iim, jjm, llm, lon, lat, lev, ttm, itaus, date0, &
       dt, fid_out)

    ! This procedure opens an input file. There is no test of the
    ! content of the file against the input from the model.

    ! The user should check that the dimensions of lon lat and lev are
    ! correct when passed to flinopen. This can be done after the call
    ! when iim and jjm have been retrieved from the netCDF file. In
    ! Fortran 90 this problem will be solved with an internal assign
    ! IF iim, jjm, llm or ttm are parameters in the calling program it
    ! will create a segmentation fault

    USE calendar, ONLY: ymds2ju, ioconf_calendar
    USE errioipsl, ONLY: histerr
    use flinfindcood_m, only: flinfindcood
    use flininfo_m, only: nbfile_max, NCFILEOPEN, NCIDS, NCNBVA
    USE netcdf, ONLY: nf90_get_att, nf90_get_var, nf90_global, &
         nf90_inquire_variable

    INTEGER, intent(in):: iim ! size in the x direction in the file (longitude)
    INTEGER, intent(in):: jjm ! size in the y direction

    INTEGER, intent(in):: llm ! number of levels
    ! (llm = 0 means no axis to be expected)

    INTEGER, intent(in):: ttm ! size of time axis
    real, intent(out):: lon(iim, jjm) ! longitude
    real, intent(out):: lat(iim, jjm) ! latitude
    real, intent(out):: lev(llm)
    INTEGER, intent(out):: itaus(ttm) ! time steps within this file
    REAL, intent(out):: date0 ! Julian date at which itau = 0
    REAL, intent(out):: dt ! length of the time steps of the data

    INTEGER, intent(in):: fid_out
    ! (file ID which is later used to read the data)

    ! Variables local to the procedure:
    INTEGER iret, vid, fid, nbdim, i
    INTEGER gdtt_id, old_id, iv, gdtmaf_id
    CHARACTER(LEN=250) name
    CHARACTER(LEN=80) units, my_calendar
    INTEGER year, month, day
    REAL r_year, r_month, r_day
    INTEGER year0, month0, day0, hours0, minutes0, seci
    REAL sec, sec0
    CHARACTER strc
    REAL, DIMENSION(:), ALLOCATABLE:: vec_tmp

    !---------------------------------------------------------------------

    IF ((fid_out < 1) .OR. (fid_out > nbfile_max)) THEN
       ! Either the fid_out has not been initialized (0 or very large)
       ! then we have to open anyway. Else we only need to open the file
       ! if it has not been opened before.
       print *, "Call flinfo before flinopen"
       stop 1
    end IF
    IF (.NOT. ncfileopen(fid_out)) THEN
       print *, "Call flinfo before flinopen"
       stop 1
    end IF

    ! The user has already opened the file
    ! and we trust that he knows the dimensions

    fid = ncids(fid_out)

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
       iret = NF90_GET_VAR (fid, vid, lat, start=(/ 1, 1 /), &
            count=(/ iim, jjm /))
    ELSE
       ALLOCATE(vec_tmp(jjm))
       iret = NF90_GET_VAR (fid, vid, vec_tmp, start=(/ 1 /), count=(/ jjm /))
       DO i=1, iim
          lat(i, :) = vec_tmp(:)
       ENDDO
       DEALLOCATE(vec_tmp)
    ENDIF

    IF (llm > 0) THEN
       CALL flinfindcood (fid_out, 'lev', vid, nbdim)
       IF (nbdim == 1) THEN
          iret = NF90_GET_VAR (fid, vid, lev, start=(/ 1 /), count=(/ llm /))
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
       IF (INDEX(my_calendar, 'XXXX') < 1) THEN
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

end module flinopen_nozoom_m
