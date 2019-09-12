module itau2date_m

  implicit none

contains

  double precision FUNCTION itau2date (itau, date0, deltat)

    ! This function transforms itau into a date. The date whith which
    ! the time axis is going to be labeled

    use calendar, only: un_jour

    INTEGER, intent(in):: itau ! current time step
    double precision, intent(in):: date0 ! Date at which itau was equal to 0
    real, intent(in):: deltat ! time step between itau s
    
    !--------------------------------------------------------------------

    itau2date = REAL(itau) * deltat / un_jour + date0

  END FUNCTION itau2date

end module itau2date_m
