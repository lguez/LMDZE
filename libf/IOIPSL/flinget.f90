MODULE flinget_m

  ! From flincom.f90, version 2.2 2006/03/07 09:21:51

  IMPLICIT NONE

  PRIVATE
  PUBLIC flinget

  INTERFACE flinget
     MODULE PROCEDURE flinget_r3d, flinget_r2d
     ! The difference between the procedures is the rank of argument "var".
  END INTERFACE

CONTAINS

  SUBROUTINE flinget_r2d(fid_in, varname, iim, jjm, llm, ttm, itau_dep, &
       itau_fin, var)

    INTEGER, intent(in):: fid_in
    CHARACTER(LEN=*), intent(in):: varname
    INTEGER, intent(in):: iim, jjm, llm, ttm, itau_dep, itau_fin
    REAL, intent(out):: var(:, :)

    ! Local:
    INTEGER :: jl, jj, ji
    REAL, DIMENSION(:), ALLOCATABLE, SAVE :: buff_tmp
    LOGICAL :: check = .FALSE.

    !---------------------------------------------------------------------

    IF (.NOT.ALLOCATED(buff_tmp)) THEN
       IF (check) WRITE(*, *) &
            "flinget_r2d : allocate buff_tmp for buff_sz = ", SIZE(var)
       ALLOCATE (buff_tmp(SIZE(var)))
    ELSE IF (SIZE(var) > SIZE(buff_tmp)) THEN
       IF (check) WRITE(*, *) &
            "flinget_r2d : re-allocate buff_tmp for buff_sz = ", SIZE(var)
       DEALLOCATE (buff_tmp)
       ALLOCATE (buff_tmp(SIZE(var)))
    ENDIF

    CALL flinget_mat(fid_in, varname, iim, jjm, llm, ttm, itau_dep, &
         itau_fin, 1, iim, 1, jjm, buff_tmp)

    jl=0
    DO jj=1, SIZE(var, 2)
       DO ji=1, SIZE(var, 1)
          jl=jl+1
          var(ji, jj) = buff_tmp(jl)
       ENDDO
    ENDDO

  END SUBROUTINE flinget_r2d

  !****************************************************************

  SUBROUTINE flinget_r3d(fid_in, varname, iim, jjm, llm, ttm, itau_dep, &
       itau_fin, var)

    INTEGER, intent(in):: fid_in
    CHARACTER(LEN=*), intent(in):: varname
    INTEGER, intent(in):: iim, jjm, llm, ttm, itau_dep, itau_fin
    REAL, intent(out):: var(:, :, :)

    ! Local:
    INTEGER :: jl, jk, jj, ji
    REAL, DIMENSION(:), ALLOCATABLE, SAVE :: buff_tmp
    LOGICAL :: check = .FALSE.

    !---------------------------------------------------------------------

    IF (.NOT.ALLOCATED(buff_tmp)) THEN
       IF (check) WRITE(*, *) &
            "flinget_r3d : allocate buff_tmp for buff_sz = ", SIZE(var)
       ALLOCATE (buff_tmp(SIZE(var)))
    ELSE IF (SIZE(var) > SIZE(buff_tmp)) THEN
       IF (check) WRITE(*, *) &
            "flinget_r3d : re-allocate buff_tmp for buff_sz = ", SIZE(var)
       DEALLOCATE (buff_tmp)
       ALLOCATE (buff_tmp(SIZE(var)))
    ENDIF

    CALL flinget_mat (fid_in, varname, iim, jjm, llm, ttm, &
         itau_dep, itau_fin, 1, iim, 1, jjm, buff_tmp)

    jl=0
    DO jk=1, SIZE(var, 3)
       DO jj=1, SIZE(var, 2)
          DO ji=1, SIZE(var, 1)
             jl=jl+1
             var(ji, jj, jk) = buff_tmp(jl)
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE flinget_r3d

  !****************************************************************

  SUBROUTINE flinget_mat(fid_in, varname, iim, jjm, llm, ttm, itau_dep, &
       itau_fin, iideb, iilen, jjdeb, jjlen, var)

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

    USE strlowercase_m,  ONLY : strlowercase
    USE errioipsl, ONLY : histerr
    USE netcdf, ONLY : nf90_byte, nf90_double, nf90_float, nf90_get_att, &
         nf90_get_var, nf90_inquire_attribute, nf90_inquire_dimension, &
         nf90_inquire_variable, nf90_inq_attname, nf90_inq_varid, nf90_int, &
         nf90_max_var_dims, nf90_noerr, nf90_short, nf90_strerror
    use flincom, only: ncids

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

    INTEGER, SAVE :: cind_vid
    INTEGER, SAVE :: cind_fid
    INTEGER, SAVE :: cind_len
    INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: cindex
    INTEGER, DIMENSION(4) :: w_sta, w_len, w_dim
    INTEGER :: iret, fid
    INTEGER :: vid, cvid, clen
    CHARACTER(LEN=70) :: str1
    CHARACTER(LEN=250) :: att_n, tmp_n
    INTEGER :: tmp_i
    REAL, SAVE :: mis_v=0.
    REAL :: tmp_r
    INTEGER :: ndims, x_typ, nb_atts
    INTEGER, DIMENSION(NF90_MAX_VAR_DIMS) :: dimids
    INTEGER :: i, i2d, cnd
    REAL, DIMENSION(:), ALLOCATABLE, SAVE :: var_tmp
    LOGICAL :: uncompress = .FALSE.
    LOGICAL :: check = .FALSE.

    !---------------------------------------------------------------------

    fid = ncids(fid_in)

    IF (check) THEN
       WRITE(*, *) &
            'flinget_mat : fid_in, fid, varname :', fid_in, fid, TRIM(varname)
       WRITE(*, *) &
            'flinget_mat : iim, jjm, llm, ttm, itau_dep, itau_fin :', &
            iim, jjm, llm, ttm, itau_dep, itau_fin
       WRITE(*, *) &
            'flinget_mat : iideb, iilen, jjdeb, jjlen :', &
            iideb, iilen, jjdeb, jjlen
    ENDIF

    uncompress = .FALSE.

    ! 1.0 We get first all the details on this variable from the file

    vid = -1
    iret = NF90_INQ_VARID (fid, varname, vid)

    IF (vid < 0 .OR. iret /= NF90_NOERR) THEN
       CALL histerr (3, 'flinget', &
            'Variable '//TRIM(varname)//' not found in file', ' ', ' ')
    ENDIF

    iret = NF90_INQUIRE_VARIABLE (fid, vid, &
         ndims=ndims, dimids=dimids, nAtts=nb_atts)
    IF (check) THEN
       WRITE(*, *) &
            'flinget_mat : fid, vid :', fid, vid
       WRITE(*, *) &
            'flinget_mat : ndims, dimids(1:ndims), nb_atts :', &
            ndims, dimids(1:ndims), nb_atts
    ENDIF

    w_dim(:) = 0
    DO i=1, ndims
       iret  = NF90_INQUIRE_DIMENSION (fid, dimids(i), len=w_dim(i))
    ENDDO
    IF (check) WRITE(*, *) &
         'flinget_mat : w_dim :', w_dim(1:ndims)

    mis_v = 0.0

    IF (nb_atts > 0) THEN
       IF (check) THEN
          WRITE(*, *) 'flinget_mat : attributes for variable :'
       ENDIF
    ENDIF
    DO i=1, nb_atts
       iret = NF90_INQ_ATTNAME (fid, vid, i, att_n)
       iret = NF90_INQUIRE_ATTRIBUTE (fid, vid, att_n, xtype=x_typ)
       CALL strlowercase (att_n)
       IF      (    (x_typ == NF90_INT).OR.(x_typ == NF90_SHORT) &
            .OR.(x_typ == NF90_BYTE) ) THEN
          iret = NF90_GET_ATT (fid, vid, att_n, tmp_i)
          IF (check) THEN
             WRITE(*, *) '   ', TRIM(att_n), ' : ', tmp_i
          ENDIF
       ELSE IF ( (x_typ == NF90_FLOAT).OR.(x_typ == NF90_DOUBLE) ) THEN
          iret = NF90_GET_ATT (fid, vid, att_n, tmp_r)
          IF (check) THEN
             WRITE(*, *) '   ', TRIM(att_n), ' : ', tmp_r
          ENDIF
          IF (index(att_n, 'missing_value') > 0) THEN
             mis_v = tmp_r
          ENDIF
       ELSE
          tmp_n = ''
          iret = NF90_GET_ATT (fid, vid, att_n, tmp_n)
          IF (check) THEN
             WRITE(*, *) '   ', TRIM(att_n), ' : ', TRIM(tmp_n)
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
       IF (check) WRITE(*, *) &
            'Dimname, iret , NF90_NOERR : ', TRIM(tmp_n), iret, NF90_NOERR
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
             iret = NF90_INQUIRE_VARIABLE (fid, cvid, xtype=x_typ, ndims=cnd)
             !-------
             IF ( cnd /= 1 .AND. x_typ /= NF90_INT) THEN
                CALL histerr (3, 'flinget', &
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
             CALL histerr (3, 'flinget', str1, &
                  'is not compressed and does not'// &
                  ' correspond to the requested size', ' ')
          ENDIF
       ELSE
          IF (w_dim(1) /= iim) THEN
             str1 = 'The longitude dimension of '//varname
             CALL histerr (3, 'flinget', str1, &
                  'in the file is not equal to the dimension', &
                  'that should be read')
          ENDIF
          IF (w_dim(2) /= jjm) THEN
             str1 = 'The latitude dimension of '//varname
             CALL histerr (3, 'flinget', str1, &
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
             CALL histerr (3, 'flinget', &
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
                CALL histerr (3, 'flinget', &
                     'You attempt to read a time slab', &
                     'but there is no time axis on this variable', varname)
             ELSE IF ((itau_fin - itau_dep) <= w_dim(i2d+2)) THEN
                w_sta(i2d+2) = itau_dep
                w_len(i2d+2) = itau_fin-itau_dep+1
             ELSE
                CALL histerr (3, 'flinget', &
                     'The time step you try to read is not', &
                     'in the file (1)', varname)
             ENDIF
          ELSE IF (ndims == i2d+2 .AND. w_dim(i2d+2) > 1) THEN
             CALL histerr (3, 'flinget', &
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
             CALL histerr (3, 'flinget', &
                  'The file contains 4 dimensions', &
                  'but only 3 are requestes for variable ', varname)
          ENDIF
          IF ((itau_fin-itau_dep) >= 0) THEN
             IF (ndims == i2d+1) THEN
                IF ((itau_fin-itau_dep) < w_dim(i2d+1) ) THEN
                   w_sta(i2d+1) = itau_dep
                   w_len(i2d+1) = itau_fin-itau_dep+1
                ELSE
                   CALL histerr (3, 'flinget', &
                        'The time step you try to read is not', &
                        'in the file (2)', varname)
                ENDIF
             ELSE
                CALL histerr (3, 'flinget', &
                     'From your input you sould have 3 dimensions', &
                     'in the file but there are 4', varname)
             ENDIF
          ELSE
             IF (ndims == i2d+1 .AND. w_dim(i2d+1) > 1) THEN
                CALL histerr (3, 'flinget', &
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

    IF (check) WRITE(*, *) &
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

    IF (check) WRITE(*, *) 'flinget_mat 3.1 : ', NF90_STRERROR (iret)

  END  SUBROUTINE flinget_mat

END MODULE flinget_m
