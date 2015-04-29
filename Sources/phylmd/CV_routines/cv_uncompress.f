module cv_uncompress_m

  implicit none

contains

  SUBROUTINE cv_uncompress(idcum, iflag, precip, cbmf, ft, fq, fu, fv, Ma, &
       qcondc, iflag1, precip1, cbmf1, ft1, fq1, fu1, fv1, Ma1, qcondc1)

    use cv_param, ONLY: nl
    USE dimphy, ONLY: klev, klon

    integer, intent(in):: idcum(:) ! (ncum)
    integer, intent(in):: iflag(klon)
    real, intent(in):: precip(klon), cbmf(klon)
    real, intent(in):: ft(klon, klev), fq(klon, klev), fu(klon, klev)
    real, intent(in):: fv(klon, klev)
    real, intent(in):: Ma(klon, klev)
    real, intent(in):: qcondc(klon, klev) !cld

    integer, intent(inout):: iflag1(klon)
    real, intent(inout):: precip1(klon), cbmf1(klon)
    real, intent(inout), dimension(klon, klev):: ft1, fq1, fu1, fv1
    real, intent(inout):: Ma1(klon, klev)
    real, intent(inout):: qcondc1(klon, klev) !cld

    ! Local:
    integer ncum, i, k

    !------------------------------------------------------------------

    ncum = size(idcum)

    do i=1, ncum
       precip1(idcum(i))=precip(i)
       cbmf1(idcum(i))=cbmf(i)
       iflag1(idcum(i))=iflag(i)
    end do

    do k=1, nl
       do i=1, ncum
          ft1(idcum(i), k)=ft(i, k)
          fq1(idcum(i), k)=fq(i, k)
          fu1(idcum(i), k)=fu(i, k)
          fv1(idcum(i), k)=fv(i, k)
          Ma1(idcum(i), k)=Ma(i, k)
          qcondc1(idcum(i), k)=qcondc(i, k)
       end do
    end do

  end SUBROUTINE cv_uncompress

end module cv_uncompress_m
