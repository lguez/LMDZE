module cv3_uncompress_m

  implicit none

contains

  SUBROUTINE cv3_uncompress(nloc, len, ncum, nd, idcum, iflag, precip, &
       VPrecip, sig, w0, ft, fq, fu, fv, inb, Ma, upwd, dnwd, dnwd0, qcondc, &
       wd, cape, da, phi, mp, iflag1, precip1, VPrecip1, sig1, w01, ft1, fq1, &
       fu1, fv1, inb1, Ma1, upwd1, dnwd1, dnwd01, qcondc1, wd1, cape1, da1, &
       phi1, mp1)

    USE cv3_param_m, ONLY: nl

    ! inputs:
    integer, intent(in):: nloc, len, ncum, nd
    integer, intent(in):: idcum(nloc)
    integer, intent(in):: iflag(nloc)
    real, intent(in):: precip(nloc)
    real, intent(in):: VPrecip(nloc, nd+1)
    real, intent(in):: sig(nloc, nd), w0(nloc, nd)
    real, intent(in):: ft(nloc, nd), fq(nloc, nd), fu(nloc, nd), fv(nloc, nd)
    integer, intent(in):: inb(nloc)
    real, intent(in):: Ma(nloc, nd)
    real, intent(in):: upwd(nloc, nd), dnwd(nloc, nd), dnwd0(nloc, nd)
    real, intent(in):: qcondc(nloc, nd)
    real, intent(in):: wd(nloc), cape(nloc)
    real, intent(in):: da(nloc, nd), phi(nloc, nd, nd), mp(nloc, nd)

    ! outputs:
    integer iflag1(len)
    real precip1(len)
    real VPrecip1(len, nd+1)
    real sig1(len, nd), w01(len, nd)
    real ft1(len, nd), fq1(len, nd), fu1(len, nd), fv1(len, nd)
    integer inb1(len)
    real Ma1(len, nd)
    real upwd1(len, nd), dnwd1(len, nd), dnwd01(len, nd)
    real qcondc1(nloc, nd)
    real wd1(nloc), cape1(nloc)
    real da1(nloc, nd), phi1(nloc, nd, nd), mp1(nloc, nd)

    ! local variables:
    integer i, k, j

    !-------------------------------------------------------------------

    do  i=1, ncum
       precip1(idcum(i))=precip(i)
       iflag1(idcum(i))=iflag(i)
       wd1(idcum(i))=wd(i)
       inb1(idcum(i))=inb(i)
       cape1(idcum(i))=cape(i)
    end do

    do  k=1, nl
       do  i=1, ncum
          VPrecip1(idcum(i), k)=VPrecip(i, k)
          sig1(idcum(i), k)=sig(i, k)
          w01(idcum(i), k)=w0(i, k)
          ft1(idcum(i), k)=ft(i, k)
          fq1(idcum(i), k)=fq(i, k)
          fu1(idcum(i), k)=fu(i, k)
          fv1(idcum(i), k)=fv(i, k)
          Ma1(idcum(i), k)=Ma(i, k)
          upwd1(idcum(i), k)=upwd(i, k)
          dnwd1(idcum(i), k)=dnwd(i, k)
          dnwd01(idcum(i), k)=dnwd0(i, k)
          qcondc1(idcum(i), k)=qcondc(i, k)
          da1(idcum(i), k)=da(i, k)
          mp1(idcum(i), k)=mp(i, k)
       end do
    end do

    do  i=1, ncum
       sig1(idcum(i), nd)=sig(i, nd)
    end do

    do j=1, nd
       do k=1, nd
          do i=1, ncum
             phi1(idcum(i), k, j)=phi(i, k, j)
          end do
       end do
    end do

  end SUBROUTINE cv3_uncompress

end module cv3_uncompress_m
