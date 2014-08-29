module cv_feed_m

  implicit none

contains

  SUBROUTINE cv_feed(len, nd, t, q, qs, p, hm, gz, nk, icb, icbmax, iflag, &
       tnk, qnk, gznk, plcl)

    use cv_param

    ! Purpose: CONVECTIVE FEED

    ! inputs:
    integer, intent(in):: len, nd
    real, intent(in):: t(len, nd)
    real, intent(in):: q(len, nd), qs(len, nd), p(len, nd)
    real hm(len, nd), gz(len, nd)

    ! outputs:
    integer iflag(len)
    integer, intent(out):: nk(len), icb(len), icbmax
    real tnk(len), qnk(len), gznk(len), plcl(len)

    ! local variables:
    integer i, k
    integer ihmin(len)
    real work(len)
    real pnk(len), qsnk(len), rh(len), chi(len)

    !-------------------------------------------------------------------
    ! --- Find level of minimum moist static energy
    ! --- If level of minimum moist static energy coincides with
    ! --- or is lower than minimum allowable parcel origin level,
    ! --- set iflag to 6.
    !-------------------------------------------------------------------

    do i=1, len
       work(i)=1.0e12
       ihmin(i)=nl
    end do
    do k=2, nlp
       do i=1, len
          if ((hm(i, k) < work(i)).and. &
               (hm(i, k) < hm(i, k-1)))then
             work(i)=hm(i, k)
             ihmin(i)=k
          endif
       end do
    end do
    do i=1, len
       ihmin(i) = min(ihmin(i), nlm)
       if (ihmin(i) <= minorig) iflag(i)=6
    end do

    !-------------------------------------------------------------------
    ! --- Find that model level below the level of minimum moist static
    ! --- energy that has the maximum value of moist static energy
    !-------------------------------------------------------------------

    do i=1, len
       work(i)=hm(i, minorig)
       nk(i)=minorig
    end do
    do k=minorig+1, nl
       do i=1, len
          if ((hm(i, k) > work(i)).and.(k <= ihmin(i)))then
             work(i)=hm(i, k)
             nk(i)=k
          endif
       end do
    end do
    !-------------------------------------------------------------------
    ! --- Check whether parcel level temperature and specific humidity
    ! --- are reasonable
    !-------------------------------------------------------------------
    do i=1, len
       if (((t(i, nk(i)) < 250.0).or. &
            (q(i, nk(i)) <= 0.0).or. &
            (p(i, ihmin(i)) < 400.0)).and. &
            (iflag(i) == 0))iflag(i)=7
    end do
    !-------------------------------------------------------------------
    ! --- Calculate lifted condensation level of air at parcel origin level
    ! --- (Within 0.2% of formula of Bolton, MON. WEA. REV., 1980)
    !-------------------------------------------------------------------
    do i=1, len
       tnk(i)=t(i, nk(i))
       qnk(i)=q(i, nk(i))
       gznk(i)=gz(i, nk(i))
       pnk(i)=p(i, nk(i))
       qsnk(i)=qs(i, nk(i))

       rh(i)=qnk(i)/qsnk(i)
       rh(i)=min(1.0, rh(i))
       chi(i)=tnk(i)/(1669.0-122.0*rh(i)-tnk(i))
       plcl(i)=pnk(i)*(rh(i)**chi(i))
       if (((plcl(i) < 200.0).or.(plcl(i) >= 2000.0)) &
            .and.(iflag(i) == 0))iflag(i)=8
    end do
    !-------------------------------------------------------------------
    ! --- Calculate first level above lcl (=icb)
    !-------------------------------------------------------------------
    do i=1, len
       icb(i)=nlm
    end do

    do k=minorig, nl
       do i=1, len
          if ((k >= (nk(i)+1)).and.(p(i, k) < plcl(i))) &
               icb(i)=min(icb(i), k)
       end do
    end do

    do i=1, len
       if ((icb(i) >= nlm).and.(iflag(i) == 0))iflag(i)=9
    end do

    ! Compute icbmax.

    icbmax=2
    do i=1, len
       icbmax=max(icbmax, icb(i))
    end do

  end SUBROUTINE cv_feed

end module cv_feed_m
