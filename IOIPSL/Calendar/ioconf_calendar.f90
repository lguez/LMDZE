module ioconf_calendar_m

  ! From IOIPSL/src/calendar.f90, version 2.0 2004/04/05 14:47:47

  ! This is the calendar used to do all calculations on time. Three
  ! types of calendars are possible:

  ! - Gregorian:
  ! The normal calendar. The time origin for the julian day in this
  ! case is 24 Nov -4713.

  ! - No leap:
  ! A 365 day year without leap years. The origin for the julian days
  ! is in this case 1 Jan 0.

  ! - xxxd:
  ! Year of xxx days with months of equal length. The origin for the
  ! julian days is then also 1 Jan 0.

  ! It is difficult to go from one calendar to the other. All
  ! operations involving julian days will be wrong. This calendar will
  ! lock the length of the year as soon as possible and forbid any
  ! further modification.

  ! For the no-leap calendar, the method is still brute force. We
  ! need to find an integer series which takes care of the length of
  ! the various month. (Jan)

  implicit none

  INTEGER:: mon_len(12) = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
  CHARACTER(LEN=20):: calendar_used = "gregorian"
  REAL:: un_an = 365.2425 ! one year in days

contains

  SUBROUTINE ioconf_calendar(str)

    ! This routine allows to configure the calendar to be used.
    ! This operation is only allowed once and the first call to
    ! ymds2ju or ju2ymsd will lock the current configuration.
    ! the argument to ioconf_calendar can be any of the following:

    ! - gregorian: this is the gregorian calendar (default here)

    ! - noleap: A calendar without leap years = 365 days

    ! - xxxd: A calendar of xxx days (has to be a modulo of 12) with
    ! 12 month of equal length

    use calendar, only: lock_unan
    use strlowercase_m, only: strlowercase
    use errioipsl, only: histerr

    CHARACTER(LEN=*), INTENT(IN):: str

    ! Local:
    INTEGER leng, ipos
    CHARACTER(LEN=10) str10
    !--------------------------------------------------------------------

    CALL strlowercase(str)

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
          ipos = INDEX(str, '_day')
          IF (ipos == 4) THEN
             READ(str(1:3), '(I3)') leng
             IF ((MOD(leng, 12) == 0).AND.(leng > 1)) THEN
                calendar_used = str
                un_an = leng
                mon_len(:) = leng
             ELSE
                CALL histerr (3, 'ioconf_calendar', &
                     'The length of the year has to be a modulo of 12', &
                     'so that it can be divided into 12 month of equal ' &
                     // 'length', str)
             ENDIF
          ELSE
             CALL histerr (3, 'ioconf_calendar', &
                  'Unrecognized input, please ceck the man pages.', str, ' ')
          ENDIF
       END SELECT
    ELSE
       WRITE(str10, '(f10.4)') un_an
       CALL histerr (2, 'ioconf_calendar', &
            'The calendar was already used or configured. You are not', &
            'allowed to change it again. '// &
            'The following length of year is used:', str10)
    ENDIF

  END SUBROUTINE ioconf_calendar

end module ioconf_calendar_m
