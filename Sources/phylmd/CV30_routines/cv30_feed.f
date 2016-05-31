module cv30_feed_m

  implicit none

contains

  SUBROUTINE cv30_feed(t1, q1, qs1, p1, ph1, gz1, icb1, iflag1, tnk1, qnk1, &
       gznk1, plcl1)

    ! Purpose: convective feed
    ! Assuming origin level of ascending parcels is minorig.

    use cv30_param_m, only: minorig, nl
    USE dimphy, ONLY: klev, klon
    use numer_rec_95, only: locate

    real, intent(in):: t1(:, :) ! (klon, klev)
    real, intent(in):: q1(:, :), qs1(:, :), p1(:, :) ! (klon, klev)
    real, intent(in):: ph1(:, :) ! (klon, klev+1)
    real, intent(in):: gz1(:, :) ! (klon, klev)

    ! outputs:

    integer, intent(out):: icb1(:) ! (klon)
    ! first level above lcl, 2 <= icb1 <= nl - 2

    integer, intent(out):: iflag1(:) ! (klon)
    real tnk1(:), qnk1(:), gznk1(:) ! (klon)
    real, intent(out):: plcl1(klon)

    ! Local:
    integer i
    real qsnk(klon), rh(klon), chi(klon)
    real, parameter:: A = 1669., B = 122.

    !--------------------------------------------------------------------

    iflag1 = 0
    plcl1 = 0.

    ! Check whether parcel level temperature and specific humidity
    ! are reasonable
    do i = 1, klon
       if (t1(i, minorig) < 250. .or. q1(i, minorig) <= 0.) iflag1(i) = 7
    end do

    ! Calculate lifted condensation level of air at parcel origin level
    ! (within 0.2 % of formula of Bolton, Mon. Wea. Rev., 1980)
    do i = 1, klon
       if (iflag1(i) == 0) then
          tnk1(i) = t1(i, minorig)
          qnk1(i) = q1(i, minorig)
          gznk1(i) = gz1(i, minorig)
          qsnk(i) = qs1(i, minorig)

          rh(i) = qnk1(i) / qsnk(i)
          chi(i) = tnk1(i) / (A - B * rh(i) - tnk1(i))
          plcl1(i) = p1(i, minorig) * (rh(i)**chi(i))
          if (plcl1(i) < 200. .or. plcl1(i) >= 2000.) iflag1(i) = 8
       endif
    end do

    ! Compute icb1:
    do i = 1, klon
       if (plcl1(i) <= ph1(i, nl - 2)) then
          ! Distinguish this case just so that icb1 = nl - 2, not nl =
          ! 3, for plcl1 exactly == ph1(i, nl - 2). Maybe not useful.
          icb1(i) = nl - 2
       else
          icb1(i) = locate(- ph1(i, 3:nl - 2), - plcl1(i), my_lbound = 3)
          ! {2 <= icb1(i) <= nl - 3}
          ! {ph1(i, icb1(i) + 1) < plcl1(i)}
          ! {plcl1(i) <= ph1(i, icb1(i)) or icb1(i) == 2}
       end if
    end do

    where (icb1 == nl - 2 .and. iflag1 == 0) iflag1 = 9

    ! {(2 <= icb1(i) <= nl - 3 and ph1(i, icb1(i) + 1) < plcl1(i) and
    ! (plcl1(i) <= ph1(i, icb1(i)) or icb1(i) == 2)) or iflag1(i) /=
    ! 0}

  end SUBROUTINE cv30_feed

end module cv30_feed_m
