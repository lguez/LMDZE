module cv30_feed_m

  implicit none

contains

  SUBROUTINE cv30_feed(len, nd, t, q, qs, p, ph, gz, nk, icb, icbmax, iflag, &
       tnk, qnk, gznk, plcl)

    ! Purpose: convective feed

    ! Main differences with cv_feed:
    ! - ph added in input
    ! - here, nk(i) = minorig
    ! - icb defined differently (plcl compared with ph instead of p)

    use cv30_param_m, only: minorig, nl, nlm

    integer, intent(in):: len, nd
    real, intent(in):: t(len, nd)
    real, intent(in):: q(len, nd), qs(len, nd), p(len, nd)
    real, intent(in):: ph(len, nd+1)
    real, intent(in):: gz(len, nd)

    ! outputs:
    integer, intent(out):: nk(len), icb(len), icbmax
    integer, intent(inout):: iflag(len)
    real tnk(len), qnk(len), gznk(len)
    real, intent(out):: plcl(len)

    ! Local:
    integer i, k
    real pnk(len), qsnk(len), rh(len), chi(len)
    real, parameter:: A = 1669., B = 122.

    !--------------------------------------------------------------------

    plcl = 0.

    ! Origin level of ascending parcels

    do i = 1, len
       nk(i) = minorig
    end do

    ! Check whether parcel level temperature and specific humidity
    ! are reasonable

    do i = 1, len
       if ((t(i, nk(i)) < 250. .or. q(i, nk(i)) <= 0.) .and. iflag(i) == 0) &
            iflag(i) = 7
    end do

    ! Calculate lifted condensation level of air at parcel origin level
    ! (within 0.2% of formula of Bolton, Mon. Wea. Rev., 1980)

    do i = 1, len
       if (iflag(i) /= 7) then
          tnk(i) = t(i, nk(i))
          qnk(i) = q(i, nk(i))
          gznk(i) = gz(i, nk(i))
          pnk(i) = p(i, nk(i))
          qsnk(i) = qs(i, nk(i))

          rh(i) = qnk(i)/qsnk(i)
          chi(i) = tnk(i)/(A-B*rh(i)-tnk(i))
          plcl(i) = pnk(i)*(rh(i)**chi(i))
          if ((plcl(i) < 200. .or. plcl(i) >= 2000.) .and. iflag(i) == 0) &
               iflag(i) = 8
       endif
    end do

    ! Calculate first level above lcl (= icb)

    do i = 1, len
       icb(i) = nlm
    end do

    ! La modification consiste \`a comparer plcl \`a ph et non \`a p:
    ! icb est d\'efini par : ph(icb) < plcl < ph(icb - 1)
    do k = 3, nl-1 ! modification pour que icb soit supérieur ou égal à 2
       do i = 1, len
          if (ph(i, k) < plcl(i)) icb(i) = min(icb(i), k)
       end do
    end do

    do i = 1, len
       if ((icb(i) == nlm).and.(iflag(i) == 0)) iflag(i) = 9
    end do

    do i = 1, len
       icb(i) = icb(i)-1 ! icb >= 2
    end do

    ! Compute icbmax

    icbmax = 2

    do i = 1, len
       if (iflag(i) < 7) icbmax = max(icbmax, icb(i))
    end do

  end SUBROUTINE cv30_feed

end module cv30_feed_m
