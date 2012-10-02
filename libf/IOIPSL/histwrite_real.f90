module histwrite_real_m

  implicit none

contains

  SUBROUTINE histwrite_real(fileid, varid, itau, nbdpt, buff_tmp, nbindex, &
       nindex, do_oper, do_write)

    ! This subroutine is internal and does the calculations and writing
    ! if needed. At a later stage it should be split into an operation
    ! and writing subroutines.

    USE histcom_var, ONLY: buffer, buff_pos, datasz_max, deltat, &
         last_opp, last_wrt, missing_val, nbopp, nb_opp, nb_wrt, ncdf_ids, &
         ncvar_ids, point, regular, scal, scsize, sopps, tax_last, tdimid, &
         topp, var_axid, zorig, zsize
    USE mathelp, ONLY: trans_buff, moycum
    USE mathop_m, ONLY: mathop
    use netcdf, only: NF90_PUT_VAR

    INTEGER, INTENT(IN):: fileid, varid, itau, nbdpt
    REAL buff_tmp(:)

    INTEGER, INTENT(IN):: nbindex
    ! number of indices provided
    ! If it is equal to the size of the full field as provided in histdef
    ! then nothing is done.

    INTEGER, INTENT(IN):: nindex(nbindex)
    ! The indices used to expand the variable (pdata) onto the full field

    LOGICAL, INTENT(IN):: do_oper, do_write

    ! Local:

    INTEGER:: tsz, ncid, ncvarid
    INTEGER:: i, iret, ipt, itax
    INTEGER:: io, nbin, nbout
    INTEGER, DIMENSION(4):: corner, edges
    INTEGER:: itime

    REAL:: rtime
    CHARACTER(LEN=7):: tmp_opp

    REAL, ALLOCATABLE, SAVE:: buff_tmp2(:)
    INTEGER, SAVE:: buff_tmp2_sz
    REAL, ALLOCATABLE, SAVE:: buffer_used(:)
    INTEGER, SAVE:: buffer_sz

    !--------------------------------------------------------------------

    ! The sizes which can be encoutered

    tsz = zsize(fileid, varid, 1) * zsize(fileid, varid, 2) &
         * zsize(fileid, varid, 3)

    ! 1.0 We allocate the memory needed to store the data between write
    ! and the temporary space needed for operations.
    ! We have to keep precedent buffer if needed

    IF (.NOT. ALLOCATED(buffer)) THEN
       ALLOCATE(buffer(buff_pos))
       buffer_sz = buff_pos
       buffer(:)=0.0
    ELSE IF (buffer_sz < buff_pos) THEN
       IF (SUM(buffer)/=0.0) THEN
          ALLOCATE(buffer_used(buffer_sz))
          buffer_used(:)=buffer(:)
          DEALLOCATE(buffer)
          ALLOCATE(buffer(buff_pos))
          buffer_sz = buff_pos
          buffer(:SIZE(buffer_used))=buffer_used
          DEALLOCATE(buffer_used)
       ELSE
          DEALLOCATE(buffer)
          ALLOCATE(buffer(buff_pos))
          buffer_sz = buff_pos
          buffer(:)=0.0
       ENDIF
    ENDIF

    ! The buffers are only deallocated when more space is needed. This
    ! reduces the umber of allocates but increases memory needs.

    IF (.NOT.ALLOCATED(buff_tmp2)) THEN
       ALLOCATE(buff_tmp2(datasz_max(fileid, varid)))
       buff_tmp2_sz = datasz_max(fileid, varid)
    ELSE IF (datasz_max(fileid, varid) > buff_tmp2_sz) THEN
       DEALLOCATE(buff_tmp2)
       ALLOCATE(buff_tmp2(datasz_max(fileid, varid)))
       buff_tmp2_sz = datasz_max(fileid, varid)
    ENDIF

    rtime = itau * deltat(fileid)
    tmp_opp = topp(fileid, varid)

    ! 3.0 Do the operations or transfer the slab of data into buff_tmp

    ! 3.1 DO the Operations only if needed

    IF (do_oper) THEN
       i = fileid
       nbout = nbdpt

       ! 3.4 We continue the sequence of operations
       ! we started in the interface routine

       DO io = 2, nbopp(i, varid), 2
          nbin = nbout
          nbout = datasz_max(i, varid)
          CALL mathop(sopps(i, varid, io), nbin, buff_tmp, missing_val, &
               nbindex, nindex, scal(i, varid, io), nbout, buff_tmp2)

          nbin = nbout
          nbout = datasz_max(i, varid)
          CALL mathop(sopps(i, varid, io+1), nbin, buff_tmp2, missing_val, &
               nbindex, nindex, scal(i, varid, io+1), nbout, buff_tmp)
       ENDDO

       ! 3.5 Zoom into the data

       CALL trans_buff(zorig(i, varid, 1), zsize(i, varid, 1), &
            zorig(i, varid, 2), zsize(i, varid, 2), zorig(i, varid, 3), &
            zsize(i, varid, 3), scsize(i, varid, 1), scsize(i, varid, 2), &
            scsize(i, varid, 3), buff_tmp, buff_tmp2_sz, buff_tmp2)

       ! 5.0 Do the operations if needed. In the case of instantaneous
       ! output we do not transfer to the buffer.

       ipt = point(fileid, varid)

       IF ((TRIM(tmp_opp) /= "inst") &
            .AND.(TRIM(tmp_opp) /= "once")) THEN
          CALL moycum(tmp_opp, tsz, buffer(ipt:), &
               buff_tmp2, nb_opp(fileid, varid))
       ENDIF

       last_opp(fileid, varid) = itau
       nb_opp(fileid, varid) = nb_opp(fileid, varid)+1

    ENDIF

    ! 6.0 Write to file if needed

    IF (do_write) THEN
       ncvarid = ncvar_ids(fileid, varid)
       ncid = ncdf_ids(fileid)

       ! 6.1 Do the operations that are needed before writting
       IF ((TRIM(tmp_opp) /= "inst") .AND. (TRIM(tmp_opp) /= "once")) THEN
          rtime = (rtime + last_wrt(fileid, varid)*deltat(fileid)) / 2.
       ENDIF

       ! 6.2 Add a value to the time axis of this variable if needed
       IF (TRIM(tmp_opp) /= "l_max" .AND. TRIM(tmp_opp) /= "l_min" &
            .AND. TRIM(tmp_opp) /= "once") THEN
          itax = var_axid(fileid, varid)
          itime = nb_wrt(fileid, varid) + 1

          IF (tax_last(fileid, itax) < itime) THEN
             iret = NF90_PUT_VAR(ncid, tdimid(fileid, itax), (/rtime/), &
                  start=(/itime/))
             tax_last(fileid, itax) = itime
          ENDIF
       ELSE
          itime=1
       ENDIF

       ! 6.3 Write the data. Only in the case of instantaneous output
       ! we do not write the buffer.

       IF (scsize(fileid, varid, 3) == 1) THEN
          IF (regular(fileid)) THEN
             corner(1:4) = (/1, 1, itime, 0/)
             edges(1:4) = (/zsize(fileid, varid, 1), &
                  zsize(fileid, varid, 2), &
                  1, 0/)
          ELSE
             corner(1:4) = (/1, itime, 0, 0/)
             edges(1:4) = (/zsize(fileid, varid, 1), 1, 0, 0/)
          ENDIF
       ELSE
          IF (regular(fileid)) THEN
             corner(1:4) = (/1, 1, 1, itime/)
             edges(1:4) = (/zsize(fileid, varid, 1), &
                  zsize(fileid, varid, 2), &
                  zsize(fileid, varid, 3), 1/)
          ELSE
             corner(1:4) = (/1, 1, itime, 0/)
             edges(1:4) = (/zsize(fileid, varid, 1), &
                  zsize(fileid, varid, 3), 1, 0/)
          ENDIF
       ENDIF

       ipt = point(fileid, varid)

       IF ((TRIM(tmp_opp) /= "inst") .AND. (TRIM(tmp_opp) /= "once")) THEN
          iret = NF90_PUT_VAR(ncid, ncvarid, buffer(ipt:), &
               start=corner(1:4), count=edges(1:4))
       ELSE
          iret = NF90_PUT_VAR(ncid, ncvarid, buff_tmp2, &
               start=corner(1:4), count=edges(1:4))
       ENDIF

       last_wrt(fileid, varid) = itau
       nb_wrt(fileid, varid) = nb_wrt(fileid, varid)+1
       nb_opp(fileid, varid) = 0
    ENDIF

  END SUBROUTINE histwrite_real

end module histwrite_real_m
