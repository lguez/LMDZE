module hgardfou_m

  IMPLICIT none

contains

  SUBROUTINE hgardfou(t, ftsol)

    ! From phylmd/hgardfou.F, v 1.1.1.1 2004/05/19 12:53:07

    ! This procedure aborts the program if the temperature gets out of range.

    USE indicesol, ONLY: nbsrf
    USE dimphy, ONLY: klev, klon
    use nr_util, only: ifirstloc

    REAL, intent(in):: t(klon, klev), ftsol(klon, nbsrf)

    ! Variables local to the procedure:

    real, parameter:: temp_min = 50., temp_max = 370. ! temperature range, in K
    INTEGER k, nsrf, jbad

    !----------------------------------------------------------

    DO k = 1, klev
       jbad = ifirstloc(t(:, k) > temp_max .or. t(:, k) < temp_min)
       if (jbad <= klon) then
          PRINT *, 'hgardfou: temperature out of range'
          print *, "t(", jbad, ", ", k, ") = ", t(jbad, k)
          stop 1
       end if
    ENDDO

    DO nsrf = 1, nbsrf
       jbad = ifirstloc(ftsol(:, nsrf) > temp_max &
            .or. ftsol(:, nsrf) < temp_min)
       if (jbad <= klon) then
          PRINT *, 'hgardfou: temperature out of range'
          print *, "ftsol(", jbad, ", ", nsrf, ") = ", ftsol(jbad, nsrf)
          stop 1
       ENDIF
    ENDDO

  END SUBROUTINE hgardfou

end module hgardfou_m
