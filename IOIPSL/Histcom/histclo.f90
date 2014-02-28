module histclo_m

  implicit none

contains

  SUBROUTINE histclo(fid)

    ! This subroutine will close the file corresponding
    ! to the argument, if the argument is present. Else it will close
    ! all opened files.

    USE errioipsl, ONLY: histerr
    USE histcom_var, ONLY: nb_files, ncdf_ids
    USE netcdf, ONLY: nf90_close, nf90_noerr

    INTEGER, INTENT (IN), OPTIONAL:: fid ! file id

    ! Variables local to the procedure:
    INTEGER ifile, ncid, iret
    INTEGER start_loop, end_loop
    CHARACTER(len=70) str70

    !---------------------------------------------------------------------

    IF (present(fid)) THEN
       start_loop = fid
       end_loop = fid
    ELSE
       start_loop = 1
       end_loop = nb_files
    END IF

    DO ifile = start_loop, end_loop
       ncid = ncdf_ids(ifile)
       iret = nf90_close(ncid)
       IF (iret/=nf90_noerr) THEN
          WRITE(str70, '("This file has already been closed:", I3)') ifile
          CALL histerr(2, 'histclo', str70, '', '')
       END IF
    END DO

  END SUBROUTINE histclo

end module histclo_m
