module itau2date_m

  implicit none

contains

  REAL FUNCTION itau2date (itau, date0, deltat)

    ! This function transforms itau into a date. The date whith which
    ! the time axis is going to be labeled

    use calendar

    ! INPUT
    !   itau: current time step
    !   date0: Date at which itau was equal to 0
    !   deltat: time step between itau s

    ! OUTPUT
    !   itau2date: Date for the given itau

    INTEGER:: itau
    REAL:: date0, deltat
    !--------------------------------------------------------------------
    itau2date = REAL(itau)*deltat/un_jour+date0

  END FUNCTION itau2date

end module itau2date_m
