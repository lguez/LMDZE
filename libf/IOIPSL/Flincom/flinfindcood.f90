module flinfindcood_m

  implicit none

contains

  SUBROUTINE flinfindcood (fid_in, axtype, vid, ndim)

    ! This subroutine explores the file in order to find
    ! the coordinate according to a number of rules

    USE errioipsl, ONLY: histerr
    use flininfo_m, only: NCIDS, NCNBVA
    USE netcdf, ONLY: nf90_get_att, nf90_inquire_dimension, &
         nf90_inquire_variable, nf90_noerr
    USE strlowercase_m, ONLY: strlowercase

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
       DO WHILE ((vid < 0).AND.(iv < ncnbva(fid_in)))
          iv = iv+1
          str1 = ''
          iret = NF90_GET_ATT (ncids(fid_in), iv, 'units', str1)
          IF (iret == NF90_NOERR) THEN
             CALL strlowercase (str1)
             IF ((INDEX(str1, TRIM(dimuni1)) == 1) &
                  .OR.(INDEX(str1, TRIM(dimuni2)) == 1) &
                  .OR.(INDEX(str1, TRIM(dimuni3)) == 1)) THEN
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
       DO WHILE ((vid < 0).AND.(iv < ncnbva(fid_in)))
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
          DO WHILE ((vid < 0).AND.(iv < ncnbva(fid_in)))
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

end module flinfindcood_m
