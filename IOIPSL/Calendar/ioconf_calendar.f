module ioconf_calendar_m

  implicit none

contains

  SUBROUTINE ioconf_calendar (str)

    ! This routine allows to configure the calendar to be used.
    ! This operation is only allowed once and the first call to
    ! ymds2ju or ju2ymsd will lock the current configuration.
    ! the argument to ioconf_calendar can be any of the following:
    !  - gregorian: This is the gregorian calendar (default here)
    !  - noleap: A calendar without leap years = 365 days
    !  - xxxd: A calendar of xxx days (has to be a modulo of 12)
    !                with 12 month of equal length

    use calendar, only: mon_len, calendar_used, lock_unan, un_an
    use strlowercase_m
    use errioipsl

    CHARACTER(LEN=*), INTENT(IN):: str

    INTEGER:: leng, ipos
    CHARACTER(LEN=10):: str10
    !--------------------------------------------------------------------

    ! 1.0 Clean up the sring !

    CALL strlowercase (str)

    IF (.NOT.lock_unan) THEN

       lock_unan=.TRUE.

       SELECT CASE(str)
       CASE('gregorian')
          calendar_used = 'gregorian'
          un_an = 365.2425
          mon_len(:)=(/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
       CASE('standard')
          calendar_used = 'gregorian'
          un_an = 365.2425
          mon_len(:)=(/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
       CASE('proleptic_gregorian')
          calendar_used = 'gregorian'
          un_an = 365.2425
          mon_len(:)=(/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
       CASE('noleap')
          calendar_used = 'noleap'
          un_an = 365.0
          mon_len(:)=(/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
       CASE('365_day')
          calendar_used = 'noleap'
          un_an = 365.0
          mon_len(:)=(/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
       CASE('365d')
          calendar_used = 'noleap'
          un_an = 365.0
          mon_len(:)=(/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
       CASE('all_leap')
          calendar_used = 'all_leap'
          un_an = 366.0
          mon_len(:)=(/31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
       CASE('366_day')
          calendar_used = 'all_leap'
          un_an = 366.0
          mon_len(:)=(/31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
       CASE('366d')
          calendar_used = 'all_leap'
          un_an = 366.0
          mon_len(:)=(/31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
       CASE DEFAULT
          ipos = INDEX(str, 'd')
          IF (ipos == 4) THEN
             READ(str(1:3), '(I3)') leng
             IF ( (MOD(leng, 12) == 0).AND.(leng > 1) ) THEN
                calendar_used = str
                un_an = leng
                mon_len(:) = leng
             ELSE
                CALL histerr (3, 'ioconf_calendar', &
                     &         'The length of the year as to be a modulo of 12', &
                     &         'so that it can be divided into 12 month of equal length', &
                     &         str)
             ENDIF
          ELSE
             CALL histerr (3, 'ioconf_calendar', &
                  &       'Unrecognized input, please ceck the man pages.', str, ' ')
          ENDIF
       END SELECT
    ELSE
       WRITE(str10, '(f10.4)') un_an
       CALL histerr (2, 'ioconf_calendar', &
            &   'The calendar was already used or configured. You are not', &
            &   'allowed to change it again. '// &
            &   'The following length of year is used:', str10)
    ENDIF

  END SUBROUTINE ioconf_calendar

end module ioconf_calendar_m
