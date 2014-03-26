module ymds2ju_m

  implicit none

contains

  SUBROUTINE ymds2ju (year, month, day, sec, julian)

    use calendar, only: un_jour
    use ymds2ju_internal_m

    INTEGER, INTENT(IN):: year, month, day
    REAL, INTENT(IN):: sec
    REAL, INTENT(OUT):: julian

    INTEGER:: julian_day
    REAL:: julian_sec

    !--------------------------------------------------------------------

    CALL ymds2ju_internal(year, month, day, sec, julian_day, julian_sec)
    julian = julian_day + julian_sec / un_jour

  END SUBROUTINE ymds2ju

end module ymds2ju_m
