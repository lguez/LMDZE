module cv30_prelim_m

  implicit none

contains

  SUBROUTINE cv30_prelim(len, nd, ndp1, t, q, p, ph, lv, cpn, tv, gz, h, hm, th)

    USE cv30_param_m, ONLY: nl
    USE cv_thermo_m, ONLY: cl, clmcpv, cpd, cpv, eps, lv0, rrd, rrv

    ! Calculate arrays of geopotential, heat capacity and static energy

    integer, intent(in):: len, nd, ndp1
    real, intent(in):: t(len, nd)
    real, intent(in):: q(len, nd)
    real, intent(in):: p(len, nd), ph(len, ndp1)

    ! outputs:
    real lv(len, nd), cpn(len, nd), tv(len, nd)
    real gz(len, nd), h(len, nd), hm(len, nd)
    real th(len, nd)

    ! Local:
    integer k, i
    real rdcp
    real tvx, tvy 
    real cpx(len, nd)

    !--------------------------------------------------------------

    do k=1, nl 
       do i=1, len
          lv(i, k)= lv0-clmcpv*(t(i, k)-273.15)
          cpn(i, k)=cpd*(1.0-q(i, k)) + cpv*q(i, k)
          cpx(i, k)=cpd*(1.0-q(i, k)) + cl*q(i, k)
          tv(i, k)=t(i, k)*(1.0 + q(i, k)/eps-q(i, k))
          rdcp=(rrd*(1.-q(i, k)) + q(i, k)*rrv)/cpn(i, k)
          th(i, k)=t(i, k)*(1000.0/p(i, k))**rdcp
       end do
    end do

    ! gz = phi at the full levels (same as p).

    do i=1, len
       gz(i, 1)=0.0
    end do

    do k=2, nl 
       do i=1, len
          tvx=t(i, k)*(1. + q(i, k)/eps-q(i, k))
          tvy=t(i, k-1)*(1. + q(i, k-1)/eps-q(i, k-1))
          gz(i, k)=gz(i, k-1) + 0.5*rrd*(tvx + tvy) &
               *(p(i, k-1)-p(i, k))/ph(i, k)
       end do
    end do

    ! h = phi + cpT (dry static energy).
    ! hm = phi + cp(T-Tbase) + Lq

    do k=1, nl 
       do i=1, len
          h(i, k)=gz(i, k) + cpn(i, k)*t(i, k)
          hm(i, k)=gz(i, k) + cpx(i, k)*(t(i, k)-t(i, 1)) + lv(i, k)*q(i, k)
       end do
    end do

  end SUBROUTINE cv30_prelim

end module cv30_prelim_m
