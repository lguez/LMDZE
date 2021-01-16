module groupeun_m

  IMPLICIT NONE

contains

  SUBROUTINE groupeun(q)

    ! From LMDZ4/libf/dyn3d/groupeun.F, version 1.1.1.1 2004/05/19
    ! 12:53:07

    USE dimensions, only: iim, llm
    USE paramet_m, only: iip1, jjp1
    USE comgeom, only: aire_2d

    REAL, intent(inout):: q(:, :, :) ! (iim + 1, :, llm)

    ! Local:
    INTEGER, PARAMETER:: ngroup=3
    REAL airen, qn
    REAL aires, qs
    INTEGER i, j, l, ig, j1, j2, i0, jd

    !-------------------------------------------------------------------

    ! Champs 3D
    jd = jjp1 - size(q, 2)
    DO l = 1, llm
       j1 = 1 + jd
       j2 = 2
       DO ig = 1, ngroup
          DO j = j1 - jd, j2 - jd
             DO i0 = 1, iim, 2**(ngroup-ig+1)
                airen = 0.
                qn = 0.
                aires = 0.
                qs = 0.
                DO i = i0, i0 + 2**(ngroup-ig+1) - 1
                   airen = airen + aire_2d(i, j)
                   aires = aires + aire_2d(i, jjp1-j+1)
                   qn = qn + q(i, j, l)
                   qs = qs + q(i, jjp1-j+1-jd, l)
                END DO
                DO i = i0, i0 + 2**(ngroup-ig+1) - 1
                   q(i, j, l) = qn*aire_2d(i, j)/airen
                   q(i, jjp1-j+1-jd, l) = qs*aire_2d(i, jjp1-j+1)/aires
                END DO
             END DO
             q(iip1, j, l) = q(1, j, l)
             q(iip1, jjp1-j+1-jd, l) = q(1, jjp1-j+1-jd, l)
          END DO
          j1 = j2 + 1
          j2 = j2 + 2**ig
       END DO
    END DO

  END SUBROUTINE groupeun

end module groupeun_m
