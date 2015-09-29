module abort_gcm_m

  implicit none

contains

  SUBROUTINE abort_gcm(modname, message)

    ! From abort_gcm.F, version 1.1.1.1 2004/05/19 12:53:05
    ! Stops the simulation, closing files and printing comments.

    USE histclo_m, only: histclo

    character(len=*), intent(in):: modname ! name of calling program
    character(len=*), intent(in):: message ! to print

    !-------------------

    print *, 'abort_gcm'
    call histclo
    print *, 'Stopping in ', modname
    print *, 'Reason: ', trim(message)
    print *, 'Houston, we have a problem.'
    STOP 1

  END SUBROUTINE abort_gcm

end module abort_gcm_m
