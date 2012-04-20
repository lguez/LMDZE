module histend_m

  implicit none

contains

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

  END SUBROUTINE histend

end module histend_m
