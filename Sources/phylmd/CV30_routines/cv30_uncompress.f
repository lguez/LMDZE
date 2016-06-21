module cv30_uncompress_m

  implicit none

contains

  SUBROUTINE cv30_uncompress(idcum, iflag, precip, VPrecip, sig, w0, ft, fq, &
       fu, fv, inb, Ma, upwd, dnwd, qcondc, cape, da, phi, mp, iflag1, &
       precip1, VPrecip1, sig1, w01, ft1, fq1, fu1, fv1, inb1, Ma1, upwd1, &
       dnwd1, qcondc1, cape1, da1, phi1, mp1)

    ! UNCOMPRESS THE FIELDS

    USE cv30_param_m, ONLY: nl
    use dimphy, only: klon, klev

    integer, intent(in):: idcum(:) ! (ncum)
    integer, intent(in):: iflag(:) ! (ncum)
    real, intent(in):: precip(klon)
    real, intent(in):: VPrecip(klon, klev+1)
    real, intent(in):: sig(klon, klev), w0(klon, klev)
    real, intent(in), dimension(klon, klev):: ft, fq, fu, fv
    integer, intent(in):: inb(:) ! (ncum)
    real, intent(in):: Ma(klon, klev)
    real, intent(in):: upwd(klon, klev), dnwd(klon, klev)
    real, intent(in):: qcondc(klon, klev)
    real, intent(in):: cape(klon)
    real, intent(in):: da(:, :) ! (klon, klev)
    real, intent(in):: phi(:, :, :) ! (klon, klev, klev)
    real, intent(in):: mp(:, :) ! (ncum, nl)

    ! outputs:
    integer, intent(out):: iflag1(:) ! (klon)
    real precip1(klon)
    real VPrecip1(klon, klev+1)
    real sig1(klon, klev), w01(klon, klev)
    real ft1(klon, klev), fq1(klon, klev), fu1(klon, klev), fv1(klon, klev)
    integer, intent(inout):: inb1(klon)
    real Ma1(klon, klev)
    real upwd1(klon, klev), dnwd1(klon, klev)
    real qcondc1(klon, klev)
    real cape1(klon)
    real, intent(inout):: da1(klon, klev), phi1(klon, klev, klev)
    real, intent(inout):: mp1(klon, klev)

    ! Local:
    integer ncum, i, k, j

    !-------------------------------------------------------------------

    ncum = size(idcum)
    iflag1 = 42 ! for non convective points

    do  i = 1, ncum
       precip1(idcum(i)) = precip(i)
       iflag1(idcum(i)) = iflag(i)
       inb1(idcum(i)) = inb(i)
       cape1(idcum(i)) = cape(i)
    end do

    do  k = 1, nl
       do  i = 1, ncum
          VPrecip1(idcum(i), k) = VPrecip(i, k)
          sig1(idcum(i), k) = sig(i, k)
          w01(idcum(i), k) = w0(i, k)
          ft1(idcum(i), k) = ft(i, k)
          fq1(idcum(i), k) = fq(i, k)
          fu1(idcum(i), k) = fu(i, k)
          fv1(idcum(i), k) = fv(i, k)
          Ma1(idcum(i), k) = Ma(i, k)
          upwd1(idcum(i), k) = upwd(i, k)
          dnwd1(idcum(i), k) = dnwd(i, k)
          qcondc1(idcum(i), k) = qcondc(i, k)
          da1(idcum(i), k) = da(i, k)
          mp1(idcum(i), k) = mp(i, k)
       end do
    end do

    do  i = 1, ncum
       sig1(idcum(i), klev) = sig(i, klev)
    end do

    do j = 1, klev
       do k = 1, klev
          do i = 1, ncum
             phi1(idcum(i), k, j) = phi(i, k, j)
          end do
       end do
    end do

  end SUBROUTINE cv30_uncompress

end module cv30_uncompress_m
