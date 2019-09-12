module isittime_m

  implicit none

contains

  SUBROUTINE isittime(itau, date0, dt, freq, last_action, last_check, do_action)

    ! This subroutine checks the time has come for a given action.
    ! This is computed from the current time-step(itau).
    ! Thus we need to have the time delta (dt), the frequency
    ! of the action (freq) and the last time it was done
    ! (last_action in units of itau).
    ! In order to extrapolate when will be the next check we need
    ! the time step of the last call (last_check).

    ! The test is done on the following condition:
    ! the distance from the current time to the time for the next
    ! action is smaller than the one from the next expected
    ! check to the next action.
    ! When the test is done on the time steps simplifactions make
    ! it more difficult to read in the code.
    ! For the real time case it is easier to understand !

    use calendar
    use itau2date_m
    use ju2ymds_m, only: ju2ymds
    use ymds2ju_m, only: ymds2ju

    INTEGER, INTENT(IN):: itau
    double precision, INTENT(IN):: date0
    REAL, INTENT(IN):: dt, freq
    INTEGER, INTENT(IN):: last_action, last_check
    LOGICAL, INTENT(OUT):: do_action

    REAL:: dt_action, dt_check
    double precision date_last_act, date_now, date_mpf, date_mp1
    double precision date_next_check, date_next_act
    INTEGER:: year, month, monthp1, day, next_check_itau, next_act_itau
    INTEGER:: yearp, dayp
    REAL:: sec, secp
    LOGICAL:: check = .FALSE.
    !--------------------------------------------------------------------
    IF (check) THEN
       WRITE(*, *) &
            &    "isittime 1.0 ", itau, date0, dt, freq, last_action, last_check
    ENDIF

    IF (last_check >= 0) THEN
       dt_action = (itau-last_action)*dt
       dt_check = (itau-last_check)*dt
       next_check_itau = itau+(itau-last_check)

       !- We are dealing with frequencies in seconds and thus operation
       !- can be done on the time steps.

       IF (freq > 0) THEN
          IF (ABS(dt_action-freq) <= ABS(dt_action+dt_check-freq)) THEN
             do_action = .TRUE.
          ELSE
             do_action = .FALSE.
          ENDIF

          !--- Here we deal with frequencies in month and work on julian days.

       ELSE
          date_now = itau2date (itau, date0, dt)
          date_last_act = itau2date (last_action, date0, dt)
          CALL ju2ymds (date_last_act, year, month, day, sec)
          monthp1 = month - freq
          yearp = year

          !--- Here we compute what logically should be the next month

          IF (month >= 13) THEN
             yearp = year+1
             monthp1 = monthp1-12
          ENDIF
          CALL ymds2ju (year, monthp1, day, sec, date_mpf)

          !--- But it could be that because of a shorter month or a bad
          !--- starting date that we end up further than we should be.
          !--- Thus we compute the first day of the next month.
          !--- We can not be beyond this date and if we are close
          !--- then we will take it as it is better.

          monthp1 = month+ABS(freq)
          yearp=year
          IF (monthp1 >= 13) THEN
             yearp = year+1
             monthp1 = monthp1 -12
          ENDIF
          dayp = 1
          secp = 0.0
          CALL ymds2ju (yearp, monthp1, dayp, secp, date_mp1)

          !--- If date_mp1 is smaller than date_mpf or only less than 4 days
          !--- larger then we take it. This needed to ensure that short month
          !--- like February do not mess up the thing !

          IF (date_mp1-date_mpf < 4.) THEN
             date_next_act = date_mp1
          ELSE
             date_next_act = date_mpf
          ENDIF
          date_next_check = itau2date (next_check_itau, date0, dt)

          !--- Transform the dates into time-steps for the needed precisions.

          next_act_itau = &
               &      last_action+INT((date_next_act-date_last_act)*(un_jour/dt))

          IF (   ABS(itau-next_act_itau) &
               &        <= ABS( next_check_itau-next_act_itau)) THEN
             do_action = .TRUE.
             IF (check) THEN
                WRITE(*, *) &
                     &         'ACT-TIME: itau, next_act_itau, next_check_itau: ', &
                     &         itau, next_act_itau, next_check_itau
                CALL ju2ymds(date_now, year, month, day, sec)
                WRITE(*, *) 'ACT-TIME: y, m, d, s: ', year, month, day, sec
                WRITE(*, *) &
                     &         'ACT-TIME: date_mp1, date_mpf: ', date_mp1, date_mpf
             ENDIF
          ELSE
             do_action = .FALSE.
          ENDIF
       ENDIF

       IF (check) THEN
          WRITE(*, *) "isittime 2.0 ", &
               &     date_next_check, date_next_act, ABS(dt_action-freq), &
               &     ABS(dt_action+dt_check-freq), dt_action, dt_check, &
               &     next_check_itau, do_action
       ENDIF
    ELSE
       do_action=.FALSE.
    ENDIF

  END SUBROUTINE isittime

end module isittime_m
