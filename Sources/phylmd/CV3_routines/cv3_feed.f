module cv3_feed_m

  implicit none

contains

  SUBROUTINE cv3_feed(len, nd, t, q, qs, p, ph, hm, gz, nk, icb, icbmax, &
       iflag, tnk, qnk, gznk, plcl)

    ! Purpose: CONVECTIVE FEED

    ! Main differences with cv_feed:
    ! - ph added in input
    ! - here, nk(i)=minorig
    ! - icb defined differently (plcl compared with ph instead of p)

    ! Main differences with convect3:
    ! - we do not compute dplcldt and dplcldr of CLIFT anymore
    ! - values iflag different (but tests identical)
    ! - A, B explicitely defined (!)

    use cv3_param_m

    ! inputs:
    integer, intent(in):: len, nd
    real, intent(in):: t(len, nd)
    real, intent(in):: q(len, nd), qs(len, nd), p(len, nd)
    real hm(len, nd), gz(len, nd)
    real, intent(in):: ph(len, nd+1)

    ! outputs:
    integer iflag(len)
    integer, intent(out):: nk(len), icb(len), icbmax
    real tnk(len), qnk(len), gznk(len), plcl(len)

    ! local variables:
    integer i, k
    integer ihmin(len)
    real work(len)
    real pnk(len), qsnk(len), rh(len), chi(len)
    real A, B ! convect3

    !--------------------------------------------------------------------

    plcl=0.0

    ! --- Origin level of ascending parcels for convect3:

    do i=1, len
       nk(i)=minorig
    end do

    ! --- Check whether parcel level temperature and specific humidity
    ! --- are reasonable

    do i=1, len
       if ((t(i, nk(i)) < 250. .or. q(i, nk(i)) <= 0.) .and. iflag(i) == 0) &
            iflag(i)=7
    end do

    ! --- Calculate lifted condensation level of air at parcel origin level
    ! --- (Within 0.2% of formula of Bolton, MON. WEA. REV., 1980)

    A = 1669.0 ! convect3
    B = 122.0 ! convect3

    do i=1, len
       if (iflag(i).ne.7) then
          tnk(i)=t(i, nk(i))
          qnk(i)=q(i, nk(i))
          gznk(i)=gz(i, nk(i))
          pnk(i)=p(i, nk(i))
          qsnk(i)=qs(i, nk(i))

          rh(i)=qnk(i)/qsnk(i)
          chi(i)=tnk(i)/(A-B*rh(i)-tnk(i)) ! convect3
          plcl(i)=pnk(i)*(rh(i)**chi(i))
          if ((plcl(i) < 200. .or. plcl(i) >= 2000.) .and. iflag(i) == 0) &
               iflag(i) = 8
       endif
    end do

    ! --- Calculate first level above lcl (=icb)

    do i=1, len
       icb(i)=nlm
    end do

    ! la modification consiste a comparer plcl a ph et non a p:
    ! icb est defini par : ph(icb) < plcl < ph(icb - 1)
    do k=3, nl-1 ! modification pour que icb soit supérieur ou égal à 2
       do i=1, len
          if(ph(i, k) < plcl(i)) icb(i) = min(icb(i), k)
       end do
    end do

    do i=1, len
       if((icb(i) == nlm).and.(iflag(i) == 0))iflag(i)=9
    end do

    do i=1, len
       icb(i) = icb(i)-1 ! icb sup ou egal a 2
    end do

    ! Compute icbmax.

    icbmax=2
    do i=1, len
       if (iflag(i) < 7) icbmax=max(icbmax, icb(i)) ! sb Jun7th02
    end do

  end SUBROUTINE cv3_feed

end module cv3_feed_m
