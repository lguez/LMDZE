MODULE ioget_calendar_m

  use calendar, only: lock_unan

  IMPLICIT NONE

  PRIVATE
  PUBLIC ioget_calendar

  INTERFACE ioget_calendar
     MODULE PROCEDURE ioget_calendar_real1, ioget_calendar_real2, &
          ioget_calendar_str
  END INTERFACE

CONTAINS

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
    use calendar, only: calendar_used

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
    use calendar, only: un_an

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
    use calendar, only: un_an, un_jour

    REAL,INTENT(OUT) :: long_an,long_jour
    !---------------------------------------------------------------------
    lock_unan = .TRUE.

    long_an = un_an
    long_jour = un_jour
    !----------------------------------
  END SUBROUTINE ioget_calendar_real2

END MODULE ioget_calendar_m
