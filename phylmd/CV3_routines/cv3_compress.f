module cv3_compress_m

  implicit none

contains

  SUBROUTINE cv3_compress(len, nloc, ncum, nd, ntra, iflag1, nk1, icb1, &
       icbs1, plcl1, tnk1, qnk1, gznk1, pbase1, buoybase1, t1, q1, qs1, u1, &
       v1, gz1, th1, tra1, h1, lv1, cpn1, p1, ph1, tv1, tp1, tvp1, clw1, &
       sig1, w01, iflag, nk, icb, icbs, plcl, tnk, qnk, gznk, pbase, &
       buoybase, t, q, qs, u, v, gz, th, tra, h, lv, cpn, p, ph, tv, tp, tvp, &
       clw, sig, w0)

    use cv3_param_m


    ! inputs:
    integer, intent(in):: len, ncum, nd, ntra, nloc
    integer iflag1(len), nk1(len), icb1(len), icbs1(len)
    real plcl1(len), tnk1(len), qnk1(len), gznk1(len)
    real pbase1(len), buoybase1(len)
    real, intent(in):: t1(len, nd)
    real, intent(in):: q1(len, nd), qs1(len, nd), u1(len, nd), v1(len, nd)
    real gz1(len, nd), h1(len, nd), lv1(len, nd), cpn1(len, nd)
    real p1(len, nd), ph1(len, nd+1), tv1(len, nd), tp1(len, nd)
    real tvp1(len, nd), clw1(len, nd)
    real th1(len, nd)
    real sig1(len, nd), w01(len, nd)
    real, intent(in):: tra1(len, nd, ntra)

    ! outputs:
    ! en fait, on a nloc=len pour l'instant (cf cv_driver)
    integer iflag(nloc), nk(nloc), icb(nloc), icbs(nloc)
    real plcl(nloc), tnk(nloc), qnk(nloc), gznk(nloc)
    real pbase(nloc), buoybase(nloc)
    real t(nloc, nd), q(nloc, nd), qs(nloc, nd), u(nloc, nd), v(nloc, nd)
    real gz(nloc, nd), h(nloc, nd), lv(nloc, nd), cpn(nloc, nd)
    real p(nloc, nd), ph(nloc, nd+1), tv(nloc, nd), tp(nloc, nd)
    real tvp(nloc, nd), clw(nloc, nd)
    real th(nloc, nd)
    real sig(nloc, nd), w0(nloc, nd)
    real tra(nloc, nd, ntra)

    ! local variables:
    integer i, k, nn, j


    do  k=1, nl+1
       nn=0
       do  i=1, len
          if(iflag1(i).eq.0)then
             nn=nn+1
             sig(nn, k)=sig1(i, k)
             w0(nn, k)=w01(i, k)
             t(nn, k)=t1(i, k)
             q(nn, k)=q1(i, k)
             qs(nn, k)=qs1(i, k)
             u(nn, k)=u1(i, k)
             v(nn, k)=v1(i, k)
             gz(nn, k)=gz1(i, k)
             h(nn, k)=h1(i, k)
             lv(nn, k)=lv1(i, k)
             cpn(nn, k)=cpn1(i, k)
             p(nn, k)=p1(i, k)
             ph(nn, k)=ph1(i, k)
             tv(nn, k)=tv1(i, k)
             tp(nn, k)=tp1(i, k)
             tvp(nn, k)=tvp1(i, k)
             clw(nn, k)=clw1(i, k)
             th(nn, k)=th1(i, k)
          endif
       end do
    end do

    if (nn.ne.ncum) then
       print*, 'strange! nn not equal to ncum: ', nn, ncum
       stop
    endif

    nn=0
    do  i=1, len
       if(iflag1(i).eq.0)then
          nn=nn+1
          pbase(nn)=pbase1(i)
          buoybase(nn)=buoybase1(i)
          plcl(nn)=plcl1(i)
          tnk(nn)=tnk1(i)
          qnk(nn)=qnk1(i)
          gznk(nn)=gznk1(i)
          nk(nn)=nk1(i)
          icb(nn)=icb1(i)
          icbs(nn)=icbs1(i)
          iflag(nn)=iflag1(i)
       endif
    end do

  end SUBROUTINE cv3_compress

end module cv3_compress_m
