MODULE calendar

  ! From IOIPSL/src/calendar.f90, version 2.0 2004/04/05 14:47:47

  !- This is the calendar used to do all calculations on time. Three
  !- types of calendars are possible :
  !-  - gregorian : The normal calendar. The time origin for the
  !-                julian day in this case is 24 Nov -4713
  !-  - nolap : A 365 day year without leap years.
  !-            The origin for the julian days is in this case 1 Jan 0
  !-  - xxxd : Year of xxx days with month of equal length.
  !-           The origin for the julian days is then also 1 Jan 0
  !- As one can see it is difficult to go from one calendar to the other.
  !- All operations involving julian days will be wrong.
  !- This calendar will lock as soon as possible
  !- the length of the year and  forbid any further modification.
  !-
  !- For the non leap-year calendar the method is still brute force.
  !- We need to find an Integer series which takes care of the length
  !- of the various month. (Jan)
  !-
  !-   un_jour : one day in seconds
  !-   un_an   : one year in days

  USE strlowercase_m, ONLY : strlowercase
  USE errioipsl, ONLY : histerr
  !-
  PRIVATE
  PUBLIC :: ymds2ju,ju2ymds,isittime,ioconf_calendar, &
       ioget_calendar,itau2date, ioconf_startdate
  !-
  INTERFACE ioget_calendar
     MODULE PROCEDURE &
          &    ioget_calendar_real1,ioget_calendar_real2,ioget_calendar_str
  END INTERFACE
  !-
  REAL,PARAMETER :: un_jour = 86400.0
  LOGICAL,SAVE :: lock_startdate = .FALSE.
  !-
  CHARACTER(LEN=30),SAVE :: time_stamp='XXXXXXXXXXXXXXXX'
  !-
  !- Description of calendar
  !-
  CHARACTER(LEN=20),SAVE :: calendar_used="gregorian"
  LOGICAL,SAVE :: lock_unan = .FALSE.
  REAL,SAVE :: un_an = 365.2425
  INTEGER,SAVE :: mon_len(12)=(/31,28,31,30,31,30,31,31,30,31,30,31/)
  !-
  !-
  !-
  CHARACTER(LEN=3),PARAMETER :: &
       &  cal(12) = (/'JAN','FEB','MAR','APR','MAY','JUN', &
       &              'JUL','AUG','SEP','OCT','NOV','DEC'/)
  !-
  REAL,SAVE :: start_day,start_sec

CONTAINS

  SUBROUTINE ymds2ju (year,month,day,sec,julian)

    IMPLICIT NONE

    INTEGER,INTENT(IN) :: year,month,day
    REAL,INTENT(IN)    :: sec
    REAL,INTENT(OUT) :: julian

    INTEGER :: julian_day
    REAL    :: julian_sec
    !---------------------------------------------------------------------
    CALL ymds2ju_internal (year,month,day,sec,julian_day,julian_sec)

    julian = julian_day + julian_sec / un_jour
    !---------------------
  END SUBROUTINE ymds2ju

  !===

  SUBROUTINE ymds2ju_internal (year,month,day,sec,julian_day,julian_sec)
    !---------------------------------------------------------------------
    !- Converts year, month, day and seconds into a julian day

    !- In 1968 in a letter to the editor of Communications of the ACM
    !- (CACM, volume 11, number 10, October 1968, p.657) Henry F. Fliegel
    !- and Thomas C. Van Flandern presented such an algorithm.

    !- See also : http://www.magnet.ch/serendipity/hermetic/cal_stud/jdn.htm

    !- In the case of the Gregorian calendar we have chosen to use
    !- the Lilian day numbers. This is the day counter which starts
    !- on the 15th October 1582.
    !- This is the day at which Pope Gregory XIII introduced the
    !- Gregorian calendar.
    !- Compared to the true Julian calendar, which starts some
    !- 7980 years ago, the Lilian days are smaler and are dealt with
    !- easily on 32 bit machines. With the true Julian days you can only
    !- the fraction of the day in the real part to a precision of
    !- a 1/4 of a day with 32 bits.
    !---------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: year,month,day
    REAL,INTENT(IN)    :: sec

    INTEGER,INTENT(OUT) :: julian_day
    REAL,INTENT(OUT)    :: julian_sec

    INTEGER :: jd,m,y,d,ml
    !---------------------------------------------------------------------
    lock_unan = .TRUE.

    m = month
    y = year
    d = day

    !- We deduce the calendar from the length of the year as it
    !- is faster than an INDEX on the calendar variable.

    !- Gregorian
    IF ( (un_an > 365.0).AND.(un_an < 366.0) ) THEN
       jd = (1461*(y+4800+INT(( m-14 )/12)))/4 &
            &      +(367*(m-2-12*(INT(( m-14 )/12))))/12 &
            &      -(3*((y+4900+INT((m-14)/12))/100))/4 &
            &      +d-32075
       jd = jd-2299160
       !- No leap or All leap
    ELSE IF (ABS(un_an-365.0) <= EPSILON(un_an) .OR. &
         &   ABS(un_an-366.0) <= EPSILON(un_an)) THEN
       ml = SUM(mon_len(1:m-1))
       jd = y*INT(un_an)+ml+(d-1)
       !- Calendar with regular month
    ELSE
       ml = INT(un_an)/12
       jd = y*INT(un_an)+(m-1)*ml+(d-1)
    ENDIF

    julian_day = jd
    julian_sec = sec
    !------------------------------
  END SUBROUTINE ymds2ju_internal
  !-
  !===
  !-
  SUBROUTINE ju2ymds (julian,year,month,day,sec)
    !---------------------------------------------------------------------
    IMPLICIT NONE

    REAL,INTENT(IN) :: julian

    INTEGER,INTENT(OUT) :: year,month,day
    REAL,INTENT(OUT)    :: sec

    INTEGER :: julian_day
    REAL    :: julian_sec
    !---------------------------------------------------------------------
    julian_day = INT(julian)
    julian_sec = (julian-julian_day)*un_jour

    CALL ju2ymds_internal(julian_day,julian_sec,year,month,day,sec)
    !---------------------
  END SUBROUTINE ju2ymds
  !-
  !===
  !-
  SUBROUTINE ju2ymds_internal (julian_day,julian_sec,year,month,day,sec)
    !---------------------------------------------------------------------
    !- This subroutine computes from the julian day the year,
    !- month, day and seconds

    !- In 1968 in a letter to the editor of Communications of the ACM
    !- (CACM, volume 11, number 10, October 1968, p.657) Henry F. Fliegel
    !- and Thomas C. Van Flandern presented such an algorithm.

    !- See also : http://www.magnet.ch/serendipity/hermetic/cal_stud/jdn.htm

    !- In the case of the Gregorian calendar we have chosen to use
    !- the Lilian day numbers. This is the day counter which starts
    !- on the 15th October 1582. This is the day at which Pope
    !- Gregory XIII introduced the Gregorian calendar.
    !- Compared to the true Julian calendar, which starts some 7980
    !- years ago, the Lilian days are smaler and are dealt with easily
    !- on 32 bit machines. With the true Julian days you can only the
    !- fraction of the day in the real part to a precision of a 1/4 of
    !- a day with 32 bits.
    !---------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: julian_day
    REAL,INTENT(IN)    :: julian_sec

    INTEGER,INTENT(OUT) :: year,month,day
    REAL,INTENT(OUT)    :: sec

    INTEGER :: l,n,i,jd,j,d,m,y,ml
    INTEGER :: add_day
    !---------------------------------------------------------------------
    lock_unan = .TRUE.

    jd = julian_day
    sec = julian_sec
    IF (sec > un_jour) THEN
       add_day = INT(sec/un_jour)
       sec = sec-add_day*un_jour
       jd = jd+add_day
    ENDIF

    !- Gregorian
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
       !- No leap or All leap
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
       !- others
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
    !------------------------------
  END SUBROUTINE ju2ymds_internal
  !-
  !===
  !-
  REAL FUNCTION itau2date (itau,date0,deltat)
    !---------------------------------------------------------------------
    !- This function transforms itau into a date. The date whith which
    !- the time axis is going to be labeled

    !- INPUT
    !-   itau   : current time step
    !-   date0  : Date at which itau was equal to 0
    !-   deltat : time step between itau s

    !- OUTPUT
    !-   itau2date : Date for the given itau
    !---------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER  :: itau
    REAL     :: date0,deltat
    !---------------------------------------------------------------------
    itau2date = REAL(itau)*deltat/un_jour+date0
    !---------------------
  END FUNCTION itau2date
  !-
  !===
  !-
  SUBROUTINE isittime &
       &  (itau,date0,dt,freq,last_action,last_check,do_action)
    !---------------------------------------------------------------------
    !- This subroutine checks the time has come for a given action.
    !- This is computed from the current time-step(itau).
    !- Thus we need to have the time delta (dt), the frequency
    !- of the action (freq) and the last time it was done
    !- (last_action in units of itau).
    !- In order to extrapolate when will be the next check we need
    !- the time step of the last call (last_check).

    !- The test is done on the following condition :
    !- the distance from the current time to the time for the next
    !- action is smaller than the one from the next expected
    !- check to the next action.
    !- When the test is done on the time steps simplifactions make
    !- it more difficult to read in the code.
    !- For the real time case it is easier to understand !
    !---------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: itau
    REAL,INTENT(IN)    :: dt,freq
    INTEGER,INTENT(IN) :: last_action,last_check
    REAL,INTENT(IN)    :: date0

    LOGICAL,INTENT(OUT)  :: do_action

    REAL :: dt_action,dt_check
    REAL :: date_last_act,date_next_check,date_next_act, &
         &        date_now,date_mp1,date_mpf
    INTEGER :: year,month,monthp1,day,next_check_itau,next_act_itau
    INTEGER :: yearp,dayp
    REAL :: sec,secp
    LOGICAL :: check = .FALSE.
    !---------------------------------------------------------------------
    IF (check) THEN
       WRITE(*,*) &
            &    "isittime 1.0 ",itau,date0,dt,freq,last_action,last_check
    ENDIF

    IF (last_check >= 0) THEN
       dt_action = (itau-last_action)*dt
       dt_check = (itau-last_check)*dt
       next_check_itau = itau+(itau-last_check)
   
       !-- We are dealing with frequencies in seconds and thus operation
       !-- can be done on the time steps.
   
       IF (freq > 0) THEN
          IF (ABS(dt_action-freq) <= ABS(dt_action+dt_check-freq)) THEN
             do_action = .TRUE.
          ELSE
             do_action = .FALSE.
          ENDIF
      
          !---- Here we deal with frequencies in month and work on julian days.
      
       ELSE
          date_now = itau2date (itau,date0,dt)
          date_last_act = itau2date (last_action,date0,dt)
          CALL ju2ymds (date_last_act,year,month,day,sec)
          monthp1 = month - freq
          yearp = year
      
          !---- Here we compute what logically should be the next month
      
          IF (month >= 13) THEN
             yearp = year+1
             monthp1 = monthp1-12
          ENDIF
          CALL ymds2ju (year,monthp1,day,sec,date_mpf)
      
          !---- But it could be that because of a shorter month or a bad
          !---- starting date that we end up further than we should be.
          !---- Thus we compute the first day of the next month.
          !---- We can not be beyond this date and if we are close
          !---- then we will take it as it is better.
      
          monthp1 = month+ABS(freq)
          yearp=year
          IF (monthp1 >= 13) THEN
             yearp = year+1
             monthp1 = monthp1 -12
          ENDIF
          dayp = 1
          secp = 0.0
          CALL ymds2ju (yearp,monthp1,dayp,secp,date_mp1)
      
          !---- If date_mp1 is smaller than date_mpf or only less than 4 days
          !---- larger then we take it. This needed to ensure that short month
          !---- like February do not mess up the thing !
      
          IF (date_mp1-date_mpf < 4.) THEN
             date_next_act = date_mp1
          ELSE
             date_next_act = date_mpf
          ENDIF
          date_next_check = itau2date (next_check_itau,date0,dt)
      
          !---- Transform the dates into time-steps for the needed precisions.
      
          next_act_itau = &
               &      last_action+INT((date_next_act-date_last_act)*(un_jour/dt))
          !-----
          IF (   ABS(itau-next_act_itau) &
               &        <= ABS( next_check_itau-next_act_itau)) THEN
             do_action = .TRUE.
             IF (check) THEN
                WRITE(*,*) &
                     &         'ACT-TIME : itau, next_act_itau, next_check_itau : ', &
                     &         itau,next_act_itau,next_check_itau
                CALL ju2ymds (date_now,year,month,day,sec)
                WRITE(*,*) 'ACT-TIME : y, m, d, s : ',year,month,day,sec
                WRITE(*,*) &
                     &         'ACT-TIME : date_mp1, date_mpf : ',date_mp1,date_mpf
             ENDIF
          ELSE
             do_action = .FALSE.
          ENDIF
       ENDIF
   
       IF (check) THEN
          WRITE(*,*) "isittime 2.0 ", &
               &     date_next_check,date_next_act,ABS(dt_action-freq), &
               &     ABS(dt_action+dt_check-freq),dt_action,dt_check, &
               &     next_check_itau,do_action
       ENDIF
    ELSE
       do_action=.FALSE.
    ENDIF
    !----------------------
  END SUBROUTINE isittime
  !-
  !===
  !-
  SUBROUTINE ioconf_calendar (str)
    !---------------------------------------------------------------------
    !- This routine allows to configure the calendar to be used.
    !- This operation is only allowed once and the first call to
    !- ymds2ju or ju2ymsd will lock the current configuration.
    !- the argument to ioconf_calendar can be any of the following :
    !-  - gregorian : This is the gregorian calendar (default here)
    !-  - noleap    : A calendar without leap years = 365 days
    !-  - xxxd      : A calendar of xxx days (has to be a modulo of 12)
    !-                with 12 month of equal length
    !---------------------------------------------------------------------
    IMPLICIT NONE

    CHARACTER(LEN=*),INTENT(IN) :: str

    INTEGER :: leng,ipos
    CHARACTER(LEN=10) :: str10
    !---------------------------------------------------------------------

    ! 1.0 Clean up the sring !

    CALL strlowercase (str)

    IF (.NOT.lock_unan) THEN
       !---
       lock_unan=.TRUE.
       !---
       SELECT CASE(str)
       CASE('gregorian')
          calendar_used = 'gregorian'
          un_an = 365.2425
          mon_len(:)=(/31,28,31,30,31,30,31,31,30,31,30,31/)
       CASE('standard')
          calendar_used = 'gregorian'
          un_an = 365.2425
          mon_len(:)=(/31,28,31,30,31,30,31,31,30,31,30,31/)
       CASE('proleptic_gregorian')
          calendar_used = 'gregorian'
          un_an = 365.2425
          mon_len(:)=(/31,28,31,30,31,30,31,31,30,31,30,31/)
       CASE('noleap')
          calendar_used = 'noleap'
          un_an = 365.0
          mon_len(:)=(/31,28,31,30,31,30,31,31,30,31,30,31/)
       CASE('365_day')
          calendar_used = 'noleap'
          un_an = 365.0
          mon_len(:)=(/31,28,31,30,31,30,31,31,30,31,30,31/)
       CASE('365d')
          calendar_used = 'noleap'
          un_an = 365.0
          mon_len(:)=(/31,28,31,30,31,30,31,31,30,31,30,31/)
       CASE('all_leap')
          calendar_used = 'all_leap'
          un_an = 366.0
          mon_len(:)=(/31,29,31,30,31,30,31,31,30,31,30,31/)
       CASE('366_day')
          calendar_used = 'all_leap'
          un_an = 366.0
          mon_len(:)=(/31,29,31,30,31,30,31,31,30,31,30,31/)
       CASE('366d')
          calendar_used = 'all_leap'
          un_an = 366.0
          mon_len(:)=(/31,29,31,30,31,30,31,31,30,31,30,31/)
       CASE DEFAULT
          ipos = INDEX(str,'d')
          IF (ipos == 4) THEN
             READ(str(1:3),'(I3)') leng
             IF ( (MOD(leng,12) == 0).AND.(leng > 1) ) THEN
                calendar_used = str
                un_an = leng
                mon_len(:) = leng
             ELSE
                CALL histerr (3,'ioconf_calendar', &
                     &         'The length of the year as to be a modulo of 12', &
                     &         'so that it can be divided into 12 month of equal length', &
                     &         str)
             ENDIF
          ELSE
             CALL histerr (3,'ioconf_calendar', &
                  &       'Unrecognized input, please ceck the man pages.',str,' ')
          ENDIF
       END SELECT
    ELSE
       WRITE(str10,'(f10.4)') un_an
       CALL histerr (2,'ioconf_calendar', &
            &   'The calendar was already used or configured. You are not', &
            &   'allowed to change it again. '// &
            &   'The following length of year is used :',str10)
    ENDIF
    !-----------------------------
  END SUBROUTINE ioconf_calendar
  !-
  !===
  !-
  SUBROUTINE ioget_calendar_str (str)
    !---------------------------------------------------------------------
    !- This subroutine returns the name of the calendar used here.
    !- Three options exist :
    !-  - gregorian : This is the gregorian calendar (default here)
    !-  - noleap    : A calendar without leap years = 365 days
    !-  - xxxd      : A calendar of xxx days (has to be a modulo of 12)
    !-                with 12 month of equal length

    !- This routine will lock the calendar.
    !- You do not want it to change after your inquiry.
    !---------------------------------------------------------------------
    IMPLICIT NONE

    CHARACTER(LEN=*),INTENT(OUT) :: str
    !---------------------------------------------------------------------
    lock_unan = .TRUE.

    str = calendar_used
    !--------------------------------
  END SUBROUTINE ioget_calendar_str
  !-
  !===
  !-
  SUBROUTINE ioget_calendar_real1 (long_an)
    !---------------------------------------------------------------------
    !- This subroutine returns the name of the calendar used here.
    !- Three options exist :
    !-  - gregorian : This is the gregorian calendar (default here)
    !-  - noleap    : A calendar without leap years = 365 days
    !-  - xxxd      : A calendar of xxx days (has to be a modulo of 12)
    !-                with 12 month of equal length

    !- This routine will lock the calendar.
    !- You do not want it to change after your inquiry.
    !---------------------------------------------------------------------
    IMPLICIT NONE

    REAL,INTENT(OUT) :: long_an
    !---------------------------------------------------------------------
    lock_unan = .TRUE.

    long_an = un_an
    !----------------------------------
  END SUBROUTINE ioget_calendar_real1
  !-
  !===
  !-
  SUBROUTINE ioget_calendar_real2 (long_an,long_jour)
    !---------------------------------------------------------------------
    !- This subroutine returns the name of the calendar used here.
    !- Three options exist :
    !-  - gregorian : This is the gregorian calendar (default here)
    !-  - noleap    : A calendar without leap years = 365 days
    !-  - xxxd      : A calendar of xxx days (has to be a modulo of 12)
    !-                with 12 month of equal length

    !- This routine will lock the calendar.
    !- You do not want it to change after your inquiry.
    !---------------------------------------------------------------------
    IMPLICIT NONE

    REAL,INTENT(OUT) :: long_an,long_jour
    !---------------------------------------------------------------------
    lock_unan = .TRUE.

    long_an = un_an
    long_jour = un_jour
    !----------------------------------
  END SUBROUTINE ioget_calendar_real2

END MODULE calendar
