module abort_gcm_m

  implicit none

contains

  SUBROUTINE abort_gcm(modname, message, ierr)

    ! From abort_gcm.F, version 1.1.1.1 2004/05/19 12:53:05

    ! Stops the simulation cleanly, closing files and printing various
    ! comments

    !  Input: modname = name of calling program
    !         message = stuff to print
    !         ierr    = severity of situation ( = 0 normal )

    USE IOIPSL, only: histclo

    character(len=*), intent(in):: modname
    integer, intent(in):: ierr
    character(len=*), intent(in):: message

    !-------------------

    print *, 'abort_gcm'

    call histclo
    print *, 'Stopping in ', modname
    print *, 'Reason = ', trim(message)
    print *, 'Houston, we have a problem ', ierr
    STOP 1

  END SUBROUTINE abort_gcm

end module abort_gcm_m
