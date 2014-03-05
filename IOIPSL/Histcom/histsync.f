module histsync_m

  implicit none

contains

  SUBROUTINE histsync(file)

    ! This subroutine will synchronise all
    ! (or one if defined) opened files.

    ! file: optional argument for fileid

    USE histcom_var, ONLY: nb_files, ncdf_ids
    USE netcdf, ONLY: nf90_sync

    INTEGER, INTENT(IN), OPTIONAL:: file

    INTEGER:: ifile, ncid, iret

    LOGICAL:: file_exists
    !---------------------------------------------------------------------

    ! 1.The loop on files to synchronise

    DO ifile = 1, nb_files

       IF (present(file)) THEN
          file_exists = (ifile==file)
       ELSE
          file_exists = .TRUE.
       END IF

       IF (file_exists) THEN
          ncid = ncdf_ids(ifile)
          iret = nf90_sync(ncid)
       END IF

    END DO

  END SUBROUTINE histsync

end module histsync_m
