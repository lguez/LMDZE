MODULE calendar

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

  ! As one can see it is difficult to go from one calendar to the
  ! other. All operations involving julian days will be wrong. This
  ! calendar will lock as soon as possible the length of the year and
  ! forbid any further modification.

  ! For the non leap-year calendar the method is still brute force.
  ! We need to find an integer series which takes care of the length
  ! of the various month. (Jan)

  IMPLICIT NONE

  REAL, PARAMETER:: un_jour = 86400. ! one day in seconds

  ! Description of calendar

  CHARACTER(LEN=20):: calendar_used = "gregorian"
  LOGICAL:: lock_unan = .FALSE.
  REAL:: un_an = 365.2425 ! one year in days
  INTEGER:: mon_len(12) = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

  CHARACTER(LEN=3), PARAMETER:: cal(12) = (/'JAN', 'FEB', 'MAR', 'APR', &
       'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'/)

  REAL, SAVE:: start_day, start_sec

END MODULE calendar
