MODULE ioget_calendar_m

  !- This subroutine returns the name of the calendar used here.
  !- Three options exist :
  !-  - gregorian : This is the gregorian calendar (default here)
  !-  - noleap    : A calendar without leap years = 365 days
  !-  - xxxd      : A calendar of xxx days (has to be a modulo of 12)
  !-                with 12 month of equal length

  !- This routine will lock the calendar.
  !- You do not want it to change after your inquiry.

  use calendar, only: lock_unan

  IMPLICIT NONE

  PRIVATE lock_unan

CONTAINS

  SUBROUTINE ioget_calendar_str (str)
    use ioconf_calendar_m, only: calendar_used

    CHARACTER(LEN=*),INTENT(OUT) :: str
    !---------------------------------------------------------------------
    lock_unan = .TRUE.

    str = calendar_used
  END SUBROUTINE ioget_calendar_str
  !-
  !===
  !-
  SUBROUTINE ioget_calendar_real(long_an,long_jour)
    use calendar, only: un_jour
    use ioconf_calendar_m, only: un_an

    REAL,INTENT(OUT) :: long_an
    REAL,INTENT(OUT), optional :: long_jour
    !---------------------------------------------------------------------
    lock_unan = .TRUE.

    long_an = un_an
    if (present(long_jour)) long_jour = un_jour
  END SUBROUTINE ioget_calendar_real

END MODULE ioget_calendar_m
