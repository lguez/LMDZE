MODULE flininfo_m

  ! From flincom.f90, version 2.2 2006/03/07 09:21:51

  IMPLICIT NONE

  ! This is the data we keep on each file we open:
  INTEGER, PARAMETER:: nbfile_max = 200
  INTEGER, SAVE:: ncids(nbfile_max)
  INTEGER, SAVE:: ncnbva(nbfile_max)
  LOGICAL, SAVE:: ncfileopen(nbfile_max)=.FALSE.

CONTAINS

  SUBROUTINE flininfo(filename, iim, jjm, llm, ttm, fid_out)

    ! This subroutine allows to get some information.
    ! It is usually called by "flinopen" but the user may want to call
    ! it before in order to allocate the space needed to extract the
    ! data from the file.

    USE errioipsl, ONLY: histerr
    USE netcdf, ONLY: nf90_inquire, nf90_inquire_dimension, nf90_nowrite
    use netcdf95, only: nf95_open
    USE strlowercase_m, ONLY: strlowercase

    CHARACTER(LEN=*), intent(in):: filename
    INTEGER, intent(out):: iim, jjm, llm, ttm, fid_out

    ! Variables local to the procedure:
    INTEGER, SAVE:: nbfiles = 0
    INTEGER, SAVE:: ncdims(nbfile_max, 4)
    INTEGER iret, fid, ndims, nvars, nb_atts, id_unlim
    INTEGER iv, lll
    CHARACTER(LEN=80) name
    CHARACTER(LEN=30) axname

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
            'Too many files. Please increase nbfile_max', &
            'in module flininfo_m.', ' ')
    ENDIF

    ncids(nbfiles) = fid
    ncdims(nbfiles, :) = (/ iim, jjm, llm, ttm /)
    ncnbva(nbfiles) = nvars
    ncfileopen(nbfiles) = .TRUE.
    fid_out = nbfiles

  END SUBROUTINE flininfo

END MODULE flininfo_m
