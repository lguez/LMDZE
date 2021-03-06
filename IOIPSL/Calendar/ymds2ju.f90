module ymds2ju_m

  implicit none

contains

  SUBROUTINE ymds2ju(year, month, day, sec, julian)

    ! Converts year, month, day and seconds into a julian day

    ! In 1968 in a letter to the editor of Communications of the ACM
    ! (CACM, volume 11, number 10, October 1968, p.657) Henry F. Fliegel
    ! and Thomas C. Van Flandern presented such an algorithm.

    ! See also: http://www.magnet.ch/serendipity/hermetic/cal_stud/jdn.htm

    ! In the case of the Gregorian calendar we have chosen to use
    ! the Lilian day numbers. This is the day counter which starts
    ! on the 15th October 1582.
    ! This is the day at which Pope Gregory XIII introduced the
    ! Gregorian calendar.
    ! Compared to the true Julian calendar, which starts some
    ! 7980 years ago, the Lilian days are smaler and are dealt with
    ! easily on 32 bit machines. With the true Julian days you can only
    ! the fraction of the day in the real part to a precision of
    ! a 1/4 of a day with 32 bits.

    USE calendar, ONLY: lock_unan, un_jour
    USE ioconf_calendar_m, ONLY: mon_len, un_an

    INTEGER, INTENT(IN):: year, month, day
    REAL, INTENT(IN):: sec
    double precision, INTENT(OUT):: julian

    ! Local:
    INTEGER jd, ml

    !--------------------------------------------------------------------

    lock_unan = .TRUE.

    ! We deduce the calendar from the length of the year as it
    ! is faster than an INDEX on the calendar variable.

    ! Gregorian
    IF ( (un_an > 365.0).AND.(un_an < 366.0) ) THEN
       jd = (1461*(year+4800+INT(( month-14 )/12)))/4 &
            +(367*(month-2-12*(INT(( month-14 )/12))))/12 &
            -(3*((year+4900+INT((month-14)/12))/100))/4 &
            +day-32075
       jd = jd-2299160
       ! No leap or All leap
    ELSE IF (ABS(un_an-365.0) <= EPSILON(un_an) .OR. &
         ABS(un_an-366.0) <= EPSILON(un_an)) THEN
       ml = SUM(mon_len(1:month-1))
       jd = year*INT(un_an)+ml+(day-1)
       ! Calendar with regular month
    ELSE
       ml = INT(un_an)/12
       jd = year*INT(un_an)+(month-1)*ml+(day-1)
    ENDIF

    julian = dble(jd) + sec / un_jour

  END SUBROUTINE ymds2ju

end module ymds2ju_m
