module ju2ymds_m

  implicit none

contains

  SUBROUTINE ju2ymds (julian, year, month, day, sec)

    ! This subroutine computes from the julian day the year,
    ! month, day and seconds

    ! In 1968 in a letter to the editor of Communications of the ACM
    ! (CACM, volume 11, number 10, October 1968, p.657) Henry F. Fliegel
    ! and Thomas C. Van Flandern presented such an algorithm.

    ! See also: http://www.magnet.ch/serendipity/hermetic/cal_stud/jdn.htm

    ! In the case of the Gregorian calendar we have chosen to use
    ! the Lilian day numbers. This is the day counter which starts
    ! on the 15th October 1582. This is the day at which Pope
    ! Gregory XIII introduced the Gregorian calendar.
    ! Compared to the true Julian calendar, which starts some 7980
    ! years ago, the Lilian days are smaler and are dealt with easily
    ! on 32 bit machines. With the true Julian days you can only the
    ! fraction of the day in the real part to a precision of a 1/4 of
    ! a day with 32 bits.

    use calendar, only: un_jour, lock_unan
    use ioconf_calendar_m, only: mon_len, un_an

    REAL, INTENT(IN):: julian
    INTEGER, INTENT(OUT):: year, month, day
    REAL, INTENT(OUT):: sec

    ! Local:
    INTEGER:: l, n, i, jd, j, d, m, y, ml
    INTEGER:: add_day

    !--------------------------------------------------------------------

    jd = INT(julian)
    sec = (julian - jd) * un_jour

    lock_unan = .TRUE.

    IF (sec > un_jour) THEN
       add_day = INT(sec / un_jour)
       sec = sec - add_day * un_jour
       jd = jd+add_day
    ENDIF

    ! Gregorian
    IF ( (un_an > 365.0).AND.(un_an < 366.0) ) THEN
       jd = jd+2299160

       l = jd+68569
       n = (4*l)/146097
       l = l-(146097*n+3)/4
       i = (4000*(l+1))/1461001
       l = l-(1461*i)/4+31
       j = (80*l)/2447
       d = l-(2447*j)/80
       l = j/11
       m = j+2-(12*l)
       y = 100*(n-49)+i+l
       ! No leap or All leap
    ELSE IF (ABS(un_an-365.0) <= EPSILON(un_an) .OR. &
         &   ABS(un_an-366.0) <= EPSILON(un_an) ) THEN
       y = jd/INT(un_an)
       l = jd-y*INT(un_an)
       m = 1
       ml = 0
       DO WHILE (ml+mon_len(m) <= l)
          ml = ml+mon_len(m)
          m = m+1
       ENDDO
       d = l-ml+1
       ! others
    ELSE
       ml = INT(un_an)/12
       y = jd/INT(un_an)
       l = jd-y*INT(un_an)
       m = (l/ml)+1
       d = l-(m-1)*ml+1
    ENDIF

    day = d
    month = m
    year = y

  END SUBROUTINE ju2ymds

end module ju2ymds_m
