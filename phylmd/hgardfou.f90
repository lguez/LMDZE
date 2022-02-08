module hgardfou_m

  IMPLICIT none

contains

  SUBROUTINE hgardfou(t_seri, ftsol)

    ! From phylmd/hgardfou.F, v 1.1.1.1, 2004/05/19 12:53:07

    ! This procedure aborts the program if the temperature gets out of range.

    use jumble, only: ifirstloc

    use abort_gcm_m, only: abort_gcm
    USE indicesol, ONLY: nbsrf, clnsurf
    USE dimphy, ONLY: klev, klon
    use phyetat0_m, only: rlon, rlat

    REAL, intent(in):: t_seri(:, :) ! (klon, klev)
    REAL, intent(in):: ftsol(:, :) ! (klon, nbsrf)

    ! Variables local to the procedure:

    real, parameter:: temp_min = 50., temp_max = 370. ! temperature range, in K
    INTEGER k, nsrf, jbad

    !----------------------------------------------------------

    DO k = 1, klev
       jbad = ifirstloc(t_seri(:, k) > temp_max .or. t_seri(:, k) < temp_min)
       if (jbad <= klon) then
          print *, "t_seri(", jbad, ", ", k, ") = ", t_seri(jbad, k)
          call abort_gcm('hgardfou', 'temperature out of range')
       end if
    ENDDO

    DO nsrf = 1, nbsrf
       jbad = ifirstloc(ftsol(:, nsrf) > temp_max &
            .or. ftsol(:, nsrf) < temp_min)
       if (jbad <= klon) then
          print *, "ftsol(position index =", jbad, ", sub-surface index =", &
               nsrf, ") =", ftsol(jbad, nsrf)
          print *, "sub-surface name: ", clnsurf(nsrf)
          print *, "longitude:", rlon(jbad), "degrees east"
          print *, "latitude:", rlat(jbad), "degrees north"
          call abort_gcm('hgardfou', 'temperature out of range')
       ENDIF
    ENDDO

  END SUBROUTINE hgardfou

end module hgardfou_m
