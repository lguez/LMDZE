module hgardfou_m

  IMPLICIT none

contains

  SUBROUTINE hgardfou(t_seri, ftsol)

    ! From phylmd/hgardfou.F, v 1.1.1.1, 2004/05/19 12:53:07

    ! This procedure aborts the program if the temperature gets out of range.

    USE indicesol, ONLY: nbsrf
    USE dimphy, ONLY: klev, klon
    use nr_util, only: ifirstloc

    REAL, intent(in):: t_seri(:, :) ! (klon, klev)
    REAL, intent(in):: ftsol(:, :) ! (klon, nbsrf)

    ! Variables local to the procedure:

    real, parameter:: temp_min = 50., temp_max = 370. ! temperature range, in K
    INTEGER k, nsrf, jbad

    !----------------------------------------------------------

    DO k = 1, klev
       jbad = ifirstloc(t_seri(:, k) > temp_max .or. t_seri(:, k) < temp_min)
       if (jbad <= klon) then
          PRINT *, 'hgardfou: temperature out of range'
          print *, "t_seri(", jbad, ", ", k, ") = ", t_seri(jbad, k)
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
