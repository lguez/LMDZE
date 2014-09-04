module histend_m

  implicit none

contains

  SUBROUTINE histend(fileid)

    ! This subroutine ends the declaration of variables, sets the time
    ! axes in the NetCDF file and puts it into write mode.

    USE errioipsl, ONLY: histerr
    USE histcom_var, ONLY: date0, freq_opp, freq_wrt, fullop, &
         missing_val, name, nb_tax, nb_var, ncdf_ids, ncvar_ids, regular, &
         tax_name, tdimid, tid, title, topp, unit_name, var_axid, var_zaxid, &
         xid, yid, zax_ids, zax_name
    USE ioget_calendar_m, ONLY: ioget_calendar_str
    USE ioipslmpp, ONLY: ioipslmpp_addatt
    USE ju2ymds_m, ONLY: ju2ymds
    USE netcdf, ONLY: nf90_float, nf90_unlimited
    use netcdf95, only: nf95_def_dim, nf95_def_var, nf95_put_att, nf95_enddef

    INTEGER, INTENT(IN):: fileid ! ID of the file to be worked on

    ! Local:
    INTEGER ncid, varid
    INTEGER ndim, iv, itx, ziv
    INTEGER itax
    INTEGER dims(4), dim_cnt
    INTEGER year, month, day, hours, minutes
    REAL sec
    REAL rtime0
    CHARACTER(len=20) tname, tunit
    CHARACTER(len=30) str30
    CHARACTER(len=80) ttitle
    CHARACTER(len=120) assoc
    CHARACTER(len=70) str70
    CHARACTER(len=3):: cal(12) = (/ 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', &
         'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'/)
    CHARACTER(len=7) tmp_opp

    !---------------------------------------------------------------------

    ncid = ncdf_ids(fileid)

    ! 1.0 Create the time axes

    call nf95_def_dim(ncid, 'time_counter', nf90_unlimited, tid(fileid))

    ! 1.1 Define all the time axes needed for this file

    DO itx = 1, nb_tax(fileid)
       IF (nb_tax(fileid)>1) THEN
          str30 = 't_' // tax_name(fileid, itx)
       ELSE
          str30 = 'time_counter'
       END IF
       call nf95_def_var(ncid, str30, nf90_float, tid(fileid), &
            tdimid(fileid, itx))

       rtime0 = date0(fileid)
       CALL ju2ymds(rtime0, year, month, day, sec)

       ! Catch any error induced by a change in calendar

       IF (year <  0) THEN
          year = 2000 + year
       END IF

       hours = int(sec/(60.*60.))
       minutes = int((sec-hours*60.*60.)/60.)
       sec = sec - (hours*60.*60.+minutes*60.)

       WRITE(str70, 7000) year, month, day, hours, minutes, int(sec)
       call nf95_put_att(ncid, tdimid(fileid, itx), 'units', trim(str70))

       CALL ioget_calendar_str(str30)
       call nf95_put_att(ncid, tdimid(fileid, itx), 'calendar', trim(str30))

       call nf95_put_att(ncid, tdimid(fileid, itx), 'title', 'Time')

       call nf95_put_att(ncid, tdimid(fileid, itx), 'long_name', &
            'Time axis')

       WRITE(str70, 7001) year, cal(month), day, hours, minutes, int(sec)
       call nf95_put_att(ncid, tdimid(fileid, itx), 'time_origin', &
            trim(str70))
    END DO

    ! 2.0 declare the variables

    DO iv = 1, nb_var(fileid)

       itax = var_axid(fileid, iv)

       tname = name(fileid, iv)
       tunit = unit_name(fileid, iv)
       ttitle = title(fileid, iv)

       IF (regular(fileid)) THEN
          dims(1:2) = (/ xid(fileid), yid(fileid) /)
          dim_cnt = 2
       ELSE
          dims(1) = xid(fileid)
          dim_cnt = 1
       END IF

       tmp_opp = topp(fileid, iv)
       ziv = var_zaxid(fileid, iv)

       ! 2.1 dimension of field

       IF ((trim(tmp_opp)/='never')) THEN
          IF ((trim(tmp_opp)/='once') .AND. (trim( &
               tmp_opp)/='l_max') .AND. (trim(tmp_opp)/='l_min')) THEN
             IF (ziv==-99) THEN
                ndim = dim_cnt + 1
                dims(dim_cnt+1:dim_cnt+2) = (/ tid(fileid), 0 /)
             ELSE
                ndim = dim_cnt + 2
                dims(dim_cnt+1:dim_cnt+2) = (/ zax_ids(fileid, ziv), &
                     tid(fileid) /)
             END IF
          ELSE
             IF (ziv==-99) THEN
                ndim = dim_cnt
                dims(dim_cnt+1:dim_cnt+2) = (/ 0, 0 /)
             ELSE
                ndim = dim_cnt + 1
                dims(dim_cnt+1:dim_cnt+2) = (/ zax_ids(fileid, ziv), 0 /)
             END IF
          END IF

          call nf95_def_var(ncid, trim(tname), nf90_float, dims(1:abs(ndim)), &
               varid)

          ncvar_ids(fileid, iv) = varid

          call nf95_put_att(ncid, varid, 'units', trim(tunit))
          call nf95_put_att(ncid, varid, 'missing_value', missing_val)
          call nf95_put_att(ncid, varid, 'long_name', trim(ttitle))
          call nf95_put_att(ncid, varid, 'short_name', trim(tname))
          call nf95_put_att(ncid, varid, 'online_operation', trim(fullop( &
               fileid, iv)))

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

          call nf95_put_att(ncid, varid, 'axis', trim(str30))

          assoc = 'nav_lat nav_lon'
          ziv = var_zaxid(fileid, iv)
          IF (ziv>0) THEN
             str30 = zax_name(fileid, ziv)
             assoc = trim(str30) // ' ' // trim(assoc)
          END IF

          IF (itax>0) THEN
             IF (nb_tax(fileid)>1) THEN
                str30 = 't_' // tax_name(fileid, itax)
             ELSE
                str30 = 'time_counter'
             END IF
             assoc = trim(str30) // ' ' // trim(assoc)

             call nf95_put_att(ncid, varid, 'interval_operation', &
                  real(freq_opp(fileid, iv)))
             call nf95_put_att(ncid, varid, 'interval_write', &
                  real(freq_wrt(fileid, iv)))
          END IF
          call nf95_put_att(ncid, varid, 'associate', trim(assoc))
       END IF
    END DO

    ! Add MPP attributes
    CALL ioipslmpp_addatt(ncid)

    ! 3.0 Put the netcdf file into write mode
    call nf95_enddef(ncid)

7000 FORMAT ('seconds since ', I4.4, '-', I2.2, '-', I2.2, ' ', I2.2, ':', &
          I2.2, ':', I2.2)
7001 FORMAT (' ', I4.4, '-', A3, '-', I2.2, ' ', I2.2, ':', I2.2, ':', I2.2)
    
  END SUBROUTINE histend

end module histend_m
