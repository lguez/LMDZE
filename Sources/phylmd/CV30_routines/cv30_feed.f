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
    ! first level above LCL, 2 <= icb1 <= nl - 2

    integer, intent(out):: iflag1(:) ! (klon)
    real tnk1(:), qnk1(:), gznk1(:) ! (klon)
    real, intent(out):: plcl1(:) ! (klon)

    ! Local:
    integer i
    real rh(klon)
    real, parameter:: A = 1669., B = 122.

    !--------------------------------------------------------------------

    ! Calculate lifted condensation level of air at parcel origin level
    ! (within 0.2 % of formula of Bolton, Mon. Wea. Rev., 1980)
    where (t1(:, minorig) >= 250. .and. q1(:, minorig) > 0.)
       ! Parcel level temperature and specific humidity are reasonable.
       tnk1 = t1(:, minorig)
       qnk1 = q1(:, minorig)
       gznk1 = gz1(:, minorig)

       rh = qnk1 / qs1(:, minorig)
       plcl1 = p1(:, minorig) * rh**(tnk1 / (A - B * rh - tnk1))
       iflag1 = 0
    elsewhere
       plcl1 = 0.
       iflag1 = 7
    end where

    where (iflag1 == 0 .and. (plcl1 < 200. .or. plcl1 >= 2000.)) iflag1 = 8

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
