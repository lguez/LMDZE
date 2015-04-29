module histvert_m

  implicit none

contains

  SUBROUTINE histvert(pfileid, pzaxname, pzaxtitle, pzaxunit, pzvalues, &
       pzaxid, pdirect)

    ! This subroutine defines a vertical axis and returns its id. It
    ! gives the user the possibility to define different vertical
    ! axes. For each variable defined with histdef a vertical axis can
    ! be specified by its ID.

    USE errioipsl, ONLY: histerr
    USE find_str_m, ONLY: find_str
    USE histcom_var, ONLY: nb_zax, nb_zax_max, ncdf_ids, zax_ids, &
         zax_name, zax_name_length, zax_size
    USE netcdf, ONLY: nf90_def_dim, nf90_def_var, nf90_enddef, &
         nf90_float, nf90_put_att, nf90_put_var, nf90_redef
    USE strlowercase_m, ONLY: strlowercase

    INTEGER, INTENT(IN):: pfileid 
    ! ID of the file the variable should be archived in

    CHARACTER(len=*), INTENT(IN):: pzaxname ! name of the vertical axis
    CHARACTER(len=*), INTENT(IN):: pzaxtitle ! title of the vertical axis
    CHARACTER(len=*), INTENT(IN):: pzaxunit ! units of the vertical axis

    REAL, INTENT(IN):: pzvalues(:) ! (pzsize) 
    ! coordinate values of the vertical axis

    INTEGER, INTENT(OUT):: pzaxid ! ID of the axis (not the netCDF ID)

    CHARACTER(len=*), INTENT(IN), OPTIONAL:: pdirect
    ! positive direction of the axis: up or down

    ! Local:
    INTEGER pzsize ! size of the vertical axis
    INTEGER:: pos, iv, nb, zdimid, zaxid_tmp
    CHARACTER(len=20):: str20, tab_str20(nb_zax_max)
    INTEGER:: tab_str20_length(nb_zax_max)
    CHARACTER(len=70):: str70, str71, str72
    CHARACTER(len=80):: str80
    CHARACTER(len=20):: direction
    INTEGER:: iret, leng, ncid

    !---------------------------------------------------------------------

    pzsize = size(pzvalues)

    ! 1.0 Verifications:
    ! Do we have enough space for an extra axis ?
    ! Is the name already in use ?

    ! - Direction of axis. Can we get if from the user.
    ! If not we put unknown.

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
       WRITE(str71, '("Check variable ", a, " in file", I3)') str20, &
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

end module histvert_m
