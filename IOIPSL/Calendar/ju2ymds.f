module ju2ymds_m

  implicit none

contains

  SUBROUTINE ju2ymds (julian, year, month, day, sec)

    use calendar, only: un_jour
    USE ju2ymds_internal_m, ONLY: ju2ymds_internal

    REAL, INTENT(IN):: julian
    INTEGER, INTENT(OUT):: year, month, day
    REAL, INTENT(OUT):: sec

    INTEGER:: julian_day
    REAL:: julian_sec
    !--------------------------------------------------------------------
    julian_day = INT(julian)
    julian_sec = (julian-julian_day)*un_jour

    CALL ju2ymds_internal(julian_day, julian_sec, year, month, day, sec)

  END SUBROUTINE ju2ymds

end module ju2ymds_m
