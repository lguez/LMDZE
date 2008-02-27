module hgardfou_m

  IMPLICIT none

contains

  SUBROUTINE hgardfou(t, tsol)

    ! From phylmd/hgardfou.F, v 1.1.1.1 2004/05/19 12:53:07

    ! This procedure aborts the program if the temperature gets out of range.

    use dimens_m
    use indicesol
    use dimphy
    use YOMCST

    REAL, intent(in):: t(klon, klev), tsol(klon, nbsrf)

    ! Variables local to the procedure:

    real, parameter:: temp_min = 50., temp_max = 370. ! temperature range, in K
    INTEGER i, k, nsrf
    INTEGER jadrs(klon), jbad
    LOGICAL ok

    !----------------------------------------------------------

    ok = .TRUE.
    DO k = 1, klev
       jbad = 0
       DO i = 1, klon
          IF (t(i, k) > temp_max) THEN
             jbad = jbad + 1
             jadrs(jbad) = i
          ENDIF
       ENDDO
       IF (jbad  >  0) THEN
          ok = .FALSE.
          DO i = 1, jbad
             print *, "t(", jadrs(i), ", ", k, ") = ", t(jadrs(i), k)
          ENDDO
       ENDIF
       jbad = 0
       DO i = 1, klon
          IF (t(i, k) < temp_min) THEN
             jbad = jbad + 1
             jadrs(jbad) = i
          ENDIF
       ENDDO
       IF (jbad  >  0) THEN
          ok = .FALSE.
          DO i = 1, jbad
             print *, "t(", jadrs(i), ", ", k, ") = ", t(jadrs(i), k)
          ENDDO
       ENDIF
    ENDDO

    DO nsrf = 1, nbsrf
       jbad = 0
       DO i = 1, klon
          IF (tsol(i, nsrf) > temp_max) THEN
             jbad = jbad + 1
             jadrs(jbad) = i
          ENDIF
       ENDDO
       IF (jbad  >  0) THEN
          ok = .FALSE.
          DO i = 1, jbad
              print *, "tsol(", jadrs(i), ", ", nsrf, ") = ", &
                  tsol(jadrs(i), nsrf)
          ENDDO
       ENDIF
       jbad = 0
       DO i = 1, klon
          IF (tsol(i, nsrf) < temp_min) THEN
             jbad = jbad + 1
             jadrs(jbad) = i
          ENDIF
       ENDDO
       IF (jbad  >  0) THEN
          ok = .FALSE.
          DO i = 1, jbad
              print *, "tsol(", jadrs(i), ", ", nsrf, ") = ", &
                  tsol(jadrs(i), nsrf)
          ENDDO
       ENDIF
    ENDDO

    IF (.NOT. ok) THEN
       PRINT *, 'hgardfou: temperature out of range'
       stop 1
    ENDIF

  END SUBROUTINE hgardfou

end module hgardfou_m
